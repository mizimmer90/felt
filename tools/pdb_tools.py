import copy
import math
import numpy as np
import mdtraj as md
import simtk.unit as unit
from pdbfixer import PDBFixer
from . import basics,minimizers

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

def create_mdtraj_from_pos(op_pos, top):
    t = md.Trajectory(op_pos/unit.nanometers, top)
    return t

def _update_mdtraj_residues(fixer_top,mdtraj_top):
    res_names = np.array(
        [res.name for res in list(fixer_top.residues())],dtype=str)
    res_ids = np.array(
        [res.id for res in list(fixer_top.residues())],dtype=int)
    for num in range(len(list(mdtraj_top.residues))):
        list(mdtraj_top.residues)[num].name = res_names[num]
        list(mdtraj_top.residues)[num].resSeq = res_ids[num]
    return mdtraj_top

def mdtraj_top_from_openmm(fixer_top):
    tmp_pdb_top = md.Topology.from_openmm(fixer_top)
    pdb_top = _update_mdtraj_residues(fixer_top,tmp_pdb_top)
    return pdb_top

def pdb_from_fixer(fixer):
    fixer_top = fixer.topology
    positions = fixer.positions
    pdb_top = mdtraj_top_from_openmm(fixer_top)
    pdb = basics.create_mdtraj_from_pos(positions,pdb_top)
    return pdb

def _fix_pdb(fixer, minimize=True, mdtraj_output=False):
    '''
    Input is an openMM object, and output is optionally an mdtraj object
    '''
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()
    if minimize:
        pdb = pdb_from_fixer(fixer)
        pdb = minimizers.minimize(pdb)
        if mdtraj_output:
            return pdb
        else:
            fixer.positions = pdb.xyz[0]
            return fixer
    else:
        if mdtraj_output:
            pdb = pdb_from_fixer(fixer)
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



