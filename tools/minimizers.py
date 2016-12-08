# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 13:31:31 2016

@author: gbowman
"""

from __future__ import print_function, division, absolute_import

import numpy as np
import simtk.unit as unit

#from . import basics
from . import basics

def minimize(struct, sim=None, max_iterations=1000, retry=True):
    # get OpenMM representation of positions (from frame 0)
    op_pos = struct.openmm_positions(0)
    
    # create a basic simulatoin if no sim object is provided
    if sim is None:
        sim = basics.setup_basic_sim(struct)
        
    if not retry:
        try:
            sim.context.setPositions(op_pos)
            sim.minimizeEnergy(maxIterations=max_iterations)
            state = sim.context.getState(getPositions=True, getEnergy=True)
            min_pos = state.getPositions(asNumpy=True)
            energy = state.getPotentialEnergy()/unit.kilojoule_per_mole
        except:
            print("Minimization failed.")
            return None
    else:
        while True:
            try:
                sim.context.setPositions(op_pos)
                sim.minimizeEnergy(maxIterations=max_iterations)
                state = sim.context.getState(getPositions=True, getEnergy=True)
                min_pos = state.getPositions(asNumpy=True)
                energy = state.getPotentialEnergy()/unit.kilojoule_per_mole
                if np.isnan(energy) or np.isnan(min_pos.sum()):
                    print("Minimization failed (NaN), trying again.")
                    continue
                else:
                    break
            except:
                print("Minimization failed, trying again.")
                
    min_struct = basics.create_mdtraj_from_pos(min_pos, struct.top)
    return min_struct

def simulated_annealing(struct, sim=None, max_T=350, min_T=200, T_spacing=0.1, steps_per_T=10):
    # get OpenMM representation of positions (from frame 0)
    op_pos = struct.openmm_positions(0)
    
    # create a basic simulatoin if no sim object is provided
    if sim is None:
        sim = basics.setup_basic_sim(struct)
        
    sim.context.setPositions(op_pos)
    for T in np.arange(max_T, min_T, -T_spacing, dtype=float):
        sim.integrator.setTemperature(T)
        sim.step(steps_per_T)
        
    state = sim.context.getState(getPositions=True)
    final_pos = state.getPositions(asNumpy=True)
    
    final_struct = basics.create_mdtraj_from_pos(final_pos, struct.top)
    return final_struct
    
def bowminimize(struct, spring_k=10.0, bottom_width=0.01, max_T=350, min_T=200, T_spacing=0.1, steps_per_T=10, prot_ff='amber03.xml', sol_ff='amber03_obc.xml', retry=True):
    # Greg Bowman's custom minimizer
    dt=0.001
        
    sim_flat_bottom = basics.setup_basic_sim_cartesian_restraint(struct, struct, spring_k, bottom_width=bottom_width, prot_ff=prot_ff, sol_ff=sol_ff, T=max_T, dt=dt)
    sim_flat_bottom.context.setVelocitiesToTemperature(max_T*unit.kelvin)
    anneal_struct = simulated_annealing(struct, sim=sim_flat_bottom, max_T=max_T, min_T=min_T, T_spacing=T_spacing, steps_per_T=steps_per_T)
    
    sim_unconstrained = basics.setup_basic_sim(anneal_struct, prot_ff=prot_ff, sol_ff=sol_ff)
    min_struct = minimize(anneal_struct, sim=sim_unconstrained, retry=retry)
    
    return min_struct

