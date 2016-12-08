import numpy as np
import mdtraj as md
import simtk.unit as unit
from pdbfixer import PDBFixer
from ..simulation_tools import basics,minimizers

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

def fix_pdb(filename):
    fixer = PDBFixer(filename)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()
    pdb = pdb_from_fixer(fixer)
    min_pdb = minimizers.minimize(pdb)
    return min_pdb
    
    
    
