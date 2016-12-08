import numpy as np
import mdtraj as md
import simtk.unit as unit
from pdbfixer import PDBFixer
from . import basics,minimizers

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








