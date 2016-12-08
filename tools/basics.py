# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 12:54:27 2016

@author: gbowman
"""

from __future__ import print_function, division, absolute_import

import mdtraj as md
import numpy as np
import scipy.spatial.distance
import simtk.openmm.app as app
import simtk.openmm as op
import simtk.unit as unit
    
def create_mdtraj_from_pos(op_pos, top):
    t = md.Trajectory(op_pos/unit.nanometers, top)
    return t
    
def setup_basic_sim(struct, prot_ff='amber03.xml', sol_ff='amber03_obc.xml', T=300, dt=0.002):
    # get OpenMM representation of topology    
    op_top = struct.top.to_openmm()
    
    forcefield = app.ForceField(prot_ff, sol_ff)
    system = forcefield.createSystem(op_top, nonbondedMethod=app.NoCutoff, constraints=None)
    integrator = op.LangevinIntegrator(T*unit.kelvin, 1./unit.picosecond, dt*unit.picoseconds)
    simulation = app.Simulation(op_top, system, integrator)

    return simulation
    
def _choose_closest_symetry_mate(start_struct, start_ind, ref_struct, ref_ind, possible_atm_names):
    start_atm = start_struct.top.atom(start_ind)
    ref_atm = ref_struct.top.atom(ref_ind)
    if start_atm.name in possible_atm_names and ref_atm.name in possible_atm_names:
        ref_res = ref_atm.residue
        atm1 = ref_struct.top.select("residue %s and name %s" % (ref_res.resSeq, possible_atm_names[0]))[0]
        atm2 = ref_struct.top.select("residue %s and name %s" % (ref_res.resSeq, possible_atm_names[1]))[0]
        d1 = scipy.spatial.distance.cdist([start_struct.xyz[0,start_atm.index]], [ref_struct.xyz[0,atm1]], metric='euclidean')
        d2 = scipy.spatial.distance.cdist([start_struct.xyz[0,start_atm.index]], [ref_struct.xyz[0,atm2]], metric='euclidean')
        if d1 <= d2:
            ref_ind = atm1
        else:
            ref_ind = atm2
        ref_atm = ref_struct.top.atom(ref_ind)
        print("Matched %s with %s" % (start_atm.name, ref_atm.name))
        
    return ref_ind

def setup_sim_cartesian_restraint_backbone(
        start_struct, ref_struct, spring_k, atm_inds=None,
        bottom_width=None, ignore_h=True, prot_ff='amber03.xml',
        sol_ff='amber03_obc.xml', T=300, dt=0.001):
    """This function is a modification of the 'setup_basic_sim_cartesian_restraint'.
       Simulations automatically constrain 'start_struct' backbone coordinates to 'ref_struct',
       but additionally constrain 'atm_inds' within 'start_struct'.
    """
    start_backbone_inds = start_struct.top.select('backbone')
    ref_backbone_inds = ref_struct.top.select('backbone')

    if atm_inds != None:
        atm_inds = np.setdiff1d(atm_inds,start_backbone_inds)

    if ignore_h:
        heavy_inds = start_struct.top.select_atom_indices("heavy")
        atm_inds = np.intersect1d(atm_inds, heavy_inds)

    if start_backbone_inds.shape[0] != ref_backbone_inds.shape[0]:
        print("ERROR: number of atom inds for the start and reference structure do not match")
        return None

    # get OpenMM representation of positions (from frame 0) and topologies
    op_top = start_struct.top.to_openmm()

    forcefield = app.ForceField(prot_ff, sol_ff)
    system = forcefield.createSystem(op_top, nonbondedMethod=app.NoCutoff, constraints=None)

    # add position restraints
    if bottom_width is None:
        force = op.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")

    else:
        force = op.CustomExternalForce("0.5*k*step((x-x0)^2+(y-y0)^2+(z-z0)^2-bottom_width*bottom_width)*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addGlobalParameter("bottom_width", bottom_width)
    force.addGlobalParameter("k", spring_k)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    # Constrain backbone
    for i in range(start_backbone_inds.shape[0]):
        start_ind = int(start_backbone_inds[i])
        ref_ind = int(ref_backbone_inds[i])

        start_atm = start_struct.top.atom(start_ind)
        atom_crd = unit.Quantity(value=ref_struct.xyz[0][ref_ind].tolist(), unit=unit.nanometers)
        force.addParticle(start_ind, atom_crd)
    # Constrain additional atoms on 'start_struct'
    if atm_inds != None:
        for i in range(atm_inds.shape[0]):
            start_ind = int(atm_inds[i])

            start_atm = start_struct.top.atom(start_ind)
            atom_crd = unit.Quantity(value=start_struct.xyz[0][start_ind].tolist(), unit=unit.nanometers)
            force.addParticle(start_ind, atom_crd)
    system.addForce(force)

    integrator = op.LangevinIntegrator(T*unit.kelvin, 1./unit.picosecond, dt*unit.picoseconds)
    simulation = app.Simulation(op_top, system, integrator)

    return simulation


    
def setup_basic_sim_cartesian_restraint(start_struct, ref_struct, spring_k, atm_inds=None, ref_atm_inds=None, bottom_width=None, ignore_h=True, prot_ff='amber03.xml', sol_ff='amber03_obc.xml', T=300, dt=0.001):
    # if no atom inds are supplied, use all atoms    
    if atm_inds is None:
        atm_inds = np.arange(start_struct.top.n_atoms)
    # if different atom inds aren't specified for the ref, then assume same as start
    if ref_atm_inds is None:
        ref_atm_inds = atm_inds
        
    if ignore_h:
        heavy_inds = start_struct.top.select_atom_indices("heavy")
        atm_inds = np.intersect1d(atm_inds, heavy_inds)
        ref_heavy_inds = ref_struct.top.select_atom_indices("heavy")
        ref_atm_inds = np.intersect1d(ref_atm_inds, ref_heavy_inds)
        
    if atm_inds.shape[0] != ref_atm_inds.shape[0]:
        print("ERROR: number of atom inds for the start and reference structure do not match")
        return None
        
    # get OpenMM representation of positions (from frame 0) and topologies
    op_top = start_struct.top.to_openmm()

    forcefield = app.ForceField(prot_ff, sol_ff)
    system = forcefield.createSystem(op_top, nonbondedMethod=app.NoCutoff, constraints=None)
    
    # add position restraints
    if bottom_width is None:
        force = op.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        
    else:
        force = op.CustomExternalForce("0.5*k*step((x-x0)^2+(y-y0)^2+(z-z0)^2-bottom_width*bottom_width)*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addGlobalParameter("bottom_width", bottom_width)
    force.addGlobalParameter("k", spring_k)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    for i in range(atm_inds.shape[0]):
        start_ind = int(atm_inds[i])
        ref_ind = int(ref_atm_inds[i])
        
        # deal with equivalence of different atoms
# WRONG! If residue shits sideways, code currently matches both atoms in start_struct to one atom in ref...
        # Arg NH1/NH2
        start_atm = start_struct.top.atom(start_ind)
#        if start_atm.residue.name is "ARG":
#            ref_ind = _choose_closest_symetry_mate(start_struct, start_ind, ref_struct, ref_ind, ['NH1', 'NH2'])
        # Asp OD1/OD2
#        elif start_atm.residue.name is "ASP":
#            ref_ind = _choose_closest_symetry_mate(start_struct, start_ind, ref_struct, ref_ind, ['OD1', 'OD2'])
        # Glu OE1/OE2
#        elif start_atm.residue.name is "GLU":
#            ref_ind = _choose_closest_symetry_mate(start_struct, start_ind, ref_struct, ref_ind, ['OE1', 'OE2'])
        # Phe CD1/CD2 and CE1/CE2
        # Tyr CD1/CD2 and CE1/CE2
#        elif start_atm.residue.name in ["PHE", "TYR"]:
#            ref_ind = _choose_closest_symetry_mate(start_struct, start_ind, ref_struct, ref_ind, ['CD1', 'CD2'])
#            ref_ind = _choose_closest_symetry_mate(start_struct, start_ind, ref_struct, ref_ind, ['CE1', 'CE2'])
        
#        if ignore_h and ref_struct.top.atom(ref_ind).element.symbol is 'H':
#            continue
        atom_crd = unit.Quantity(value=ref_struct.xyz[0][ref_ind].tolist(), unit=unit.nanometers)
        force.addParticle(start_ind, atom_crd)
    system.addForce(force)

    integrator = op.LangevinIntegrator(T*unit.kelvin, 1./unit.picosecond, dt*unit.picoseconds)
    simulation = app.Simulation(op_top, system, integrator)
    
    return simulation
    
def setup_basic_sim_torsion_restraint(start_struct, ref_struct, spring_k, atm_inds=None, ref_atm_inds=None, bottom_width=None, ignore_h=True, prot_ff='amber03.xml', sol_ff='amber03_obc.xml', T=300, dt=0.001):
    # if no atom inds are supplied, use all dihedrals
    use_all_dihedrals =  False    
    if atm_inds is None:
        use_all_dihedrals = True
    # if different atom inds aren't specified for the ref, then assume same as start
    elif ref_atm_inds is None:
        ref_atm_inds = atm_inds    
    elif atm_inds.shape[0] != ref_atm_inds.shape[0]:
        print("ERROR: number of atom inds for the start and reference structure do not match")
        return None
        
# CHECK BOTTOM_WIDTH IN BOUNDS
# CHECK ATM_INDS IS RIGHT SHAPE
        
    # get OpenMM representation of positions (from frame 0) and topologies
    op_top = start_struct.top.to_openmm()

    forcefield = app.ForceField(prot_ff, sol_ff)
    system = forcefield.createSystem(op_top, nonbondedMethod=app.NoCutoff, constraints=None)
    
    # add position restraints
    if bottom_width is None:
        force = op.CustomTorsionForce("0.5*k*((cos(theta)-cos(theta0))^2+(sin(theta)-sin(theta0))^2)")
    else:
        force = op.CustomTorsionForce("0.5*k*step((cos(theta)-cos(theta0))^2+(sin(theta)-sin(theta0))^2-bottom_width*bottom_width)*((cos(theta)-cos(theta0))^2+(sin(theta)-sin(theta0))^2)")
        force.addGlobalParameter("bottom_width", bottom_width)
    force.addGlobalParameter("k", spring_k)
    force.addPerTorsionParameter("theta0")
    
    if use_all_dihedrals:
        all_dihedral_methods = [md.compute_phi, md.compute_psi, md.compute_chi1, md.compute_chi2, md.compute_chi3, md.compute_chi4, md.compute_omega]
        for dihedral_method in all_dihedral_methods:
            atms, dihedrals = dihedral_method(ref_struct)
            dihedrals = dihedrals[0]
            n_dihedrals = dihedrals.shape[0]
            for i in range(n_dihedrals):
                force.addTorsion(int(atms[i][0]), int(atms[i][1]), int(atms[i][2]), int(atms[i][3]), [dihedrals[i]])
    else:
        dihedrals = md.compute_dihedrals(ref_struct, ref_atm_inds)[0]
        n_dihedrals = len(atm_inds)
        for i in range(n_dihedrals):
            atms = atm_inds[i]
            force.addTorsion(int(atms[0]), int(atms[1]), int(atms[2]), int(atms[3]), [dihedrals[i]])
    system.addForce(force)

    integrator = op.LangevinIntegrator(T*unit.kelvin, 1./unit.picosecond, dt*unit.picoseconds)
    simulation = app.Simulation(op_top, system, integrator)
    
    return simulation
    
def get_energy(struct, sim=None):
     # get OpenMM representation of positions (from frame 0)
    op_pos = struct.openmm_positions(0)
    
    # create a basic simulatoin if no sim object is provided
    if sim is None:
        sim = setup_basic_sim(struct)
        
    sim.context.setPositions(op_pos)
    state = sim.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()/unit.kilojoule_per_mole
    
    return energy
    
def get_traj_energies(traj, sim=None):
    # create a basic simulatoin if no sim object is provided
    if sim is None:
        sim = setup_basic_sim(traj[0])
        
    n_structs = len(traj)
    energies = np.inf*np.ones(n_structs)
    for i in range(n_structs):
        energies[i] = get_energy(traj[i], sim=sim)
        
    return energies
    
def get_traj_barrier(traj, sim=None):
    energies = get_traj_energies(traj, sim=sim)
    
    # error if starting energy is nan
    if np.isnan(energies[0]):
        print("ERROR: starting structure in traj has nan energy")
        return energies[0]
        
    # make all relative to starting energy
    #energies -= energies[0]
    
    # ignore nan values
    energies = energies[np.where(np.isnan(energies)==False)[0]]
    
    return energies.max()
    
