# Authors: Maxwell I. Zimmerman, Gregory R. Bowman
# Copyright (c) 2016, Washington University
# Washington University in St. Louis

##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division, absolute_import

import numpy as np
import simtk.unit as unit
from . import sim_basics

##############################################################################
# Code
##############################################################################

def minimize(struct, sim=None, max_iterations=1000,retry=True):
    # get OpenMM representation of positions (from frame 0)
    op_pos = struct.openmm_positions(0)
    # create a basic simulatoin if no sim object is provided
    if sim is None:
        sim = sim_basics.setup_basic_sim(struct)
    while True:
        try:
            sim.context.setPositions(op_pos)
            sim.minimizeEnergy(maxIterations=max_iterations)
            state = sim.context.getState(getPositions=True, getEnergy=True)
            min_pos = state.getPositions(asNumpy=True)
            energy = state.getPotentialEnergy()/unit.kilojoule_per_mole
            if np.isnan(energy) and not retry:
                print("Minimization failed (NaN).")
                break
            elif np.isnan(energy):
                print("Minimization failed (NaN), trying again.")
            else:
                break
        except:
            print("Minimization failed.")
            if not retry:
                break
    min_struct = sim_basics.create_mdtraj_from_pos(min_pos,struct.top)
    return min_struct

def simulated_annealing(
        struct, sim=None, max_T=350, min_T=200, T_spacing=0.1, steps_per_T=10):
    # get OpenMM representation of positions (from frame 0)
    op_pos = struct.openmm_positions(0)
    # create a basic simulatoin if no sim object is provided
    if sim is None:
        sim = sim_basics.setup_basic_sim(struct)
    sim.context.setPositions(op_pos)
    for T in np.arange(max_T, min_T, -T_spacing, dtype=float):
        sim.integrator.setTemperature(T)
        sim.step(steps_per_T)
    state = sim.context.getState(getPositions=True)
    final_pos = state.getPositions(asNumpy=True)
    final_struct = sim_basics.create_mdtraj_from_pos(final_pos, struct.top)
    return final_struct

def anneal_fixer_sidechains(
        fixer, spring_const=15., bottom_width=0.005, T_max=350, T_min=200,
        T_spacing=0.1, steps_per_T=15, dt=0.002, prot_ff='amber03.xml',
        sol_ff='amber03_obc.xml'):
    """Anneals a fixer while restraining backbone atoms
    """
    # Setup ff params for minimization
    pdb = sim_basics.pdb_from_fixer(fixer)
    sim_setup = sim_basics.setup_basic_sim(pdb,prot_ff=prot_ff, sol_ff=sol_ff)
    # Minimize
    min_pdb1 = minimize(pdb, sim=sim_setup)
    # Anneal
    backbone_iis = pdb.topology.select("backbone")
    sim_flat_bottom = sim_basics.setup_basic_sim_cartesian_restraint(
        min_pdb1, min_pdb1, spring_const, atm_inds=backbone_iis,
        bottom_width=bottom_width, dt=dt, prot_ff=prot_ff, sol_ff=sol_ff)
    annealed = simulated_annealing(
        min_pdb1, sim=sim_flat_bottom, max_T=T_max, min_T=T_min,
        T_spacing=T_spacing, steps_per_T=steps_per_T)
    # Minimize p2 and output
    min_pdb2 = minimize(annealed, sim=sim_setup)
    fixer.positions = min_pdb2.openmm_positions(0)
    energy = sim_basics.get_energy(min_pdb2,sim=sim_setup)
    return fixer,energy

def _partial_restrain_fixer_from_ref(
        fixer, pdb_ref, steps=1000, spring_const=15., bottom_width=0.005, dt=0.002,
        iis_struct=None, prot_ff='amber03.xml', sol_ff='amber03_obc.xml'):
    """Constrains the backbone of fixer to the coordinates of pdb_ref,
       and additionally restrains iis_struct within fixer.
       fixer is an openMM object, pdb_ref is an MDTraj object"""
    pdb = sim_basics.pdb_from_fixer(fixer)
    sim_flat_bottom = sim_basics.setup_sim_cartesian_restraint_backbone(
        pdb, pdb_ref, spring_const, bottom_width=bottom_width, dt=dt,
        atm_inds=iis_struct, prot_ff=prot_ff, sol_ff=sol_ff)
    sim_flat_bottom.context.setPositions(fixer.positions)
    sim_flat_bottom.step(int(steps))
    state = sim_flat_bottom.context.getState(getPositions=True)
    fixer.positions = state.getPositions(asNumpy=True)
    return fixer

def relax_mutation(
        fixer, pdb_ref, max_iters1=200, max_iters2=2000, md_steps=1000,
        spring_const=15., bottom_width=0.005, dt=0.002, iis_struct=None,
        prot_ff='amber03.xml', sol_ff='amber03_obc.xml'):
    """Performs a simulation on fixer by constraining backbone atoms to
       pdb_ref and additionally restraining iis_struct. This function
       is intended to be used after mutating a residue to settle sidechains
       and surrounding atoms into place. In this example, iis_struct would
       be all of the atom indices that are further than some distance from
       the mutation site."""
    pdb = sim_basics.pdb_from_fixer(fixer)
    sim_setup = sim_basics.setup_basic_sim(pdb,prot_ff=prot_ff, sol_ff=sol_ff)
    # Minimize
    min_pdb1 = minimize(pdb, sim=sim_setup, max_iterations=max_iters1)
    fixer.positions = min_pdb1.openmm_positions(0)
    # Relax structure
    relaxed_fixer = _partial_restrain_fixer_from_ref(
        fixer, pdb_ref, steps=md_steps, spring_const=spring_const,
        bottom_width=bottom_width, dt=dt, iis_struct=iis_struct,
        prot_ff=prot_ff, sol_ff=sol_ff)
    # Minimize p2 and output
    pdb = sim_basics.pdb_from_fixer(relaxed_fixer)
    min_pdb2 = minimize(pdb, sim=sim_setup, max_iterations=max_iters2)
    fixer.positions = min_pdb2.openmm_positions(0)
    energy = sim_basics.get_energy(min_pdb2,sim=sim_setup)
    return fixer,energy


