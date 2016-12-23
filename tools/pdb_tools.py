# Author: Maxwell I. Zimmerman
# Washington University in St. Louis

###########################################################
# Imports
###########################################################

import copy
import math
import numpy as np
import mdtraj as md
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

#def rotate_sidechain(pdb, 
