# Author: Maxwell I. Zimmerman
# Washington University in St. Louis

###########################################################
# Imports
###########################################################

import copy
import math
import numpy as np
import mdtraj as md
import scipy.linalg
import simtk.unit as unit
from pdbfixer import PDBFixer
from . import sim_basics, minimizers

###########################################################
# Code
###########################################################

convert_map_1letter = {
    'A' : 'ALA', 'N' : 'ASN', 'C' : 'CYS', 'E' : 'GLU', 'H' : 'HIS',
    'L' : 'LEU', 'M' : 'MET', 'P' : 'PRO', 'T' : 'THR', 'Y' : 'TYR',
    'R' : 'ARG', 'D' : 'ASP', 'Q' : 'GLN', 'G' : 'GLY', 'I' : 'ILE',
    'K' : 'LYS', 'F' : 'PHE', 'S' : 'SER', 'W' : 'TRP', 'V' : 'VAL'}

convert_map_3letter = {
    'ALA' : 'A', 'ASN' : 'N', 'CYS' : 'C', 'GLU' : 'E', 'HIS' : 'H',
    'LEU' : 'L', 'MET' : 'M', 'PRO' : 'P', 'THR' : 'T', 'TYR' : 'Y',
    'ARG' : 'R', 'ASP' : 'D', 'GLN' : 'Q', 'GLY' : 'G', 'ILE' : 'I',
    'LYS' : 'K', 'PHE' : 'F', 'SER' : 'S', 'TRP' : 'W', 'VAL' : 'V'}

protein_residues_3letter = [
    'ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR',
    'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']

protein_residues_1letter = [
    'A', 'N', 'C', 'E', 'H', 'L', 'M', 'P', 'T', 'Y',
    'R', 'D', 'Q', 'G', 'I', 'K', 'F', 'S', 'W', 'V']


def select_random_mutation(exclude=False):
    """Provides a random amino acid 3-letter code. Amino acids can be
       excluded from being selected. For Example:

    >>> select_random_mutation(exclude=['ALA','LEU'])

    """
    if exclude:
        protein_residues = np.setdiff1d(protein_residues_3letter, exclude)
    return protein_residues[np.random.randint(len(protein_residues))]

def convert_1letter_seq(seq):
    """Converts from a 1 letter amino acid to the 3 letter code"""
    seq = [convert_map_1letter[aa] for aa in seq]
    return seq

def convert_3letter_seq(seq, concat_output=False):
    """Converts a list of 3-letter amino acid seq into a single list.
       
       >>> convert_3letter_seq(['ALA','PHE','GLY'])
       >>> 'AFG'
    
    """
    new_seq = [convert_map_3letter[aa] for aa in seq]
    if concat_output:
        new_seq = "".join(new_seq)
    return new_seq

def get_sequence(fixer,res_subset=None):
    current_residues = list(fixer.topology.residues())
    full_sequence = np.array([res.name for res in current_residues])
    if res_subset is None:
        current_seq_3letter = full_sequence
    else:
        full_res_nums = np.array([int(res.id) for res in current_residues])
        current_seq_3letter = []
        for num in res_subset:
            ii = np.where(full_res_nums==num)
            current_seq_3letter.append(full_sequence[ii][0])
    return current_seq_3letter

def _fix_pdb(fixer, minimize=True, mdtraj_output=False):
    '''
    Input is an openMM object, and output is optionally an mdtraj object
    '''
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()
    if minimize:
        pdb = sim_basics.pdb_from_fixer(fixer)
        pdb = minimizers.minimize(pdb)
        if mdtraj_output:
            return pdb
        else:
            fixer.positions = pdb.xyz[0]
            return fixer
    else:
        if mdtraj_output:
            pdb = sim_basics.pdb_from_fixer(fixer)
            return pdb
        else:
            return fixer

def fix_pdb(filename): 
    fixer = PDBFixer(filename)
    pdb = _fix_pdb(fixer, minimize=True, mdtraj_output=True)
    return pdb

def _apply_mutations(fixer, change_list, max_attempts=5):
    residue_list = list(fixer.topology.residues())
    chain_id = list(fixer.topology.chains())[0].id
    success = False
    attempts = 0
    while attempts < max_attempts and not success:
        fixer_copy = copy.deepcopy(fixer)
        fixer_copy.applyMutations(change_list, chain_id)
        try:
            fixer_copy = _fix_pdb(
                fixer_copy, minimize=True, mdtraj_output=False)
            if not math.isnan(fixer_copy.positions[0][0]._value):
                success = True
            else:
                raise
        except:
            print("failed to mutate residue...\tAttempt: "+str(attempts))
        attemps += 1
    return fixer_copy, success

def rodrigues_rotation(v, k, theta, center=None):
    if center is None:
        center = np.array([0,0,0])
    v_centered = v-center
    first_term = v_centered*np.cos(theta)
    second_term = np.cross(k,v_centered)*np.sin(theta)
    third_term = k*np.dot(k,v_centered)*(1-np.cos(theta))
    rotated_v = first_term+second_term+third_term
    return rotated_v+center

def __get_iis_and_vec(pdb, res_num):
    '''This is a helper function to rotate_chi1. Given the pdb and residue
       number, it will return the indices of atoms to rotate, the normalized
       vector to rotate along, and the center-coordinate for rotation
       reference (the beta-carbon coordinate).
    '''
    # Get the atomic indices of the atoms to rotate
    top = pdb.topology
    coords = pdb.xyz[0]
    backbone_names = ['CA','CB','N','O','H','HA']
    md_search_string = "name "+" or name ".join(backbone_names)
    backbone_iis = top.select(md_search_string)
    res_iis = top.select("resSeq %d" % res_num)
    rotate_iis = np.setdiff1d(res_iis,backbone_iis)
    # Get the vector of rotation (alpha-carbon to beta-carbon)
    alpha_ii = top.select("resSeq %d and name CA" % res_num)[0]
    beta_ii = top.select("resSeq %d and name CB" % res_num)[0]
    alpha_coord = coords[alpha_ii]
    beta_coord = coords[beta_ii]
    unnormed_rot_vec = beta_coord-alpha_coord
    rotation_vec = unnormed_rot_vec/scipy.linalg.norm(unnormed_rot_vec)
    return rotate_iis,rotation_vec,beta_coord

def rotate_chi1(pdb, res_num, thetas=None):
    top = pdb.topology
    res_nums = np.array([res.resSeq for res in top.residues])
    res_names = np.array([res.name for res in top.residues])
    if not np.any(res_nums==res_num):
        raise
    res_ii = np.where(res_nums==res_num)
    res_name = res_names[res_ii]
    if res_name == 'PRO':
        raise
    if res_name == 'GLY':
        raise
    if res_name == 'ALA':
        raise
    rotate_iis,rotation_vec,center_coord = __get_iis_and_vec(pdb,res_num)
    if thetas is None:
        thetas = [math.pi*2/3., np.pi*4/3.]
    new_pdbs = [pdb]
    for theta in thetas:
        new_pdb = copy.deepcopy(pdb)
        for ii in rotate_iis:
            new_pdb.xyz[0][ii] = rodrigues_rotation(
                new_pdb.xyz[0][ii], rotation_vec, theta, center=center_coord)
        new_pdb = minimizers.minimize(new_pdb,max_iterations=100)
        new_pdbs.append(new_pdb)
    return new_pdbs
        



