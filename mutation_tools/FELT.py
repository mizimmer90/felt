# Free-Energy Landscape Tool (FELT)
# Author: Maxwell I. Zimmerman
# Washington University in St. Louis

###########################################################
# Imports
###########################################################

from ..simulation_tools import basics
import copy
import gc
import math
import mdtraj as md
from ..simulation_tools import minimizers
import numpy as np
import os
import subprocess as sp
import sys
import time
from itertools import repeat
from multiprocessing import Pool
from pdbfixer import PDBFixer
from sabres.analysis.PageRanks import PageRank
from sabres.analysis.Rankings import FAST_Rankings
from simtk.openmm.app import PDBFile

###########################################################
# Functions
###########################################################

def select_random_mutation(exclude=False):
    """Provides a random amino acid 3-letter code. Amino acids can be
       excluded from being selected. For Example:

    >>> select_random_mutation(exclude=['ALA','LEU'])

    """
    proteinResidues = [
        'ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR',
        'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']
    if exclude:
        proteinResidues = np.setdiff1d(proteinResidues,exclude)
    return proteinResidues[np.random.randint(len(proteinResidues))]

def convert_3letter_code(seq):
    """Converts a list of 3-letter amino acid codes into a single list.
       
       >>> convert_3letter_code(['ALA','PHE','GLY'])
       >>> ['A','F','G']
    
    """
    convert_map = {
        'ALA' : 'A', 'ASN' : 'N', 'CYS' : 'C', 'GLU' : 'E', 'HIS' : 'H',
        'LEU' : 'L', 'MET' : 'M', 'PRO' : 'P', 'THR' : 'T', 'TYR' : 'Y',
        'ARG' : 'R', 'ASP' : 'D', 'GLN' : 'Q', 'GLY' : 'G', 'ILE' : 'I',
        'LYS' : 'K', 'PHE' : 'F', 'SER' : 'S', 'TRP' : 'W', 'VAL' : 'V'}
    new_seq = [convert_map[aa] for aa in seq]
    return new_seq

def get_pdb_filenames(path, absolute_pathname=True):
    """Finds the filenames that are within a
    specified directory.
    """
    if absolute_pathname:
        home_directory = sp.check_output("pwd")[:-1].decode('utf-8')+"/"
        os.chdir(path)
        path_directory = sp.check_output("pwd")[:-1].decode('utf-8')+"/"
        filenames = sp.check_output(["ls"]).decode('utf-8').split()
    else:
        filenames = sp.check_output(["ls",path]).decode('utf-8').split()
    pdb_filenames = []
    for filename in filenames:
        if filename[-4:] == '.pdb':
            if absolute_pathname:
                pdb_filenames.append(path_directory+filename)
            else:
                pdb_filenames.append(filename)
    if absolute_pathname:
        os.chdir(home_directory)
    return np.array(pdb_filenames,dtype=str)

def get_residue_list(pdb,allowed_residues=False):
    fixer = PDBFixer(filename=pdb)
    res_nums = np.array(
        [int(res.id) for res in list(fixer.topology.residues())],dtype=int)
    res_names = np.array(
        [str(res.name) for res in list(fixer.topology.residues())],dtype=str)
    if np.any(allowed_residues):
        iis = np.where(
            np.sum(
                [res_nums==res_num for res_num in allowed_residues], axis=0))
    else:
        iis = np.where(res_nums>=0)
    res_list = np.array(
        list(zip(res_nums[iis],res_names[iis])),
        dtype=[('number','int'),('name',np.str_,3)])
    return res_list

def check_topologies(pdbs):
    pdb_list = md.load(pdbs[0])
    for pdb in pdbs[1:]:
        try:
            pdb_list = pdb_list.join(md.load(pdb))
        except:
            print("\nTopology of "+str(pdb.split("/")[-1])+" does not match!")
            sys.exit()
    return pdb_list

def load_fixers(pdbs,output_directory="",mut_num=None):
    fixers = []
    for pdb in pdbs:
        if mut_num == None:
            directory = output_directory+""
        else:
            directory = output_directory+"mutation_num"+str(mut_num)+"/"
        fixer = PDBFixer(filename=directory+pdb)
        fixers.append(fixer)
    return fixers

def check_prev_seqs(prev_seqs,new_seq_1d):
    all_seqs = copy.copy(prev_seqs)
    all_seqs.append(new_seq_1d)
    if len(all_seqs) == len(np.unique(all_seqs)):
        seq_previously_used = False
    else:
        seq_previously_used = True
        print("Sequence was used already")
    return seq_previously_used

def generate_random_mutations(
        residue_numbers, current_seq, prev_seqs=False, verbose=True):
    if prev_seqs:
        seq_previously_used = True
        error_check = 0
        while seq_previously_used:
            rand_ii = np.random.randint(0,len(residue_numbers))
            rand_mutation = select_random_mutation(exclude=current_seq[rand_ii])
            new_seq = current_seq.copy()
            new_seq[rand_ii] = rand_mutation
            new_seq_1d = "".join(convert_3letter_code(new_seq))
            seq_previously_used = check_prev_seqs(prev_seqs,new_seq_1d)
            error_check += 1
            if error_check > 100:
                print("Something went wrong...")
                sys.exit()
    else:
        rand_ii = np.random.randint(0,len(residue_numbers))
        rand_mutation = select_random_mutation(exclude=current_seq[rand_ii])
        new_seq = current_seq.copy()
        new_seq[rand_ii] = rand_mutation
        new_seq_1d = "".join(convert_3letter_code(new_seq))
    change_list = current_seq[rand_ii]+\
        "-"+str(residue_numbers[rand_ii])+"-"+rand_mutation
    if verbose:
        print(
            "  Mutating residue "+
            str(current_seq[rand_ii])+"-"+
            str(residue_numbers[rand_ii])+" to "+rand_mutation)
    return [change_list], new_seq, new_seq_1d

def apply_mutation(fixer, change_list, max_iterations=5):
    residue_list = list(fixer.topology.residues())
    chain_id = list(fixer.topology.chains())[0].id
    fixer_copy = copy.deepcopy(fixer)
    print("        -mutating residue")
    fixer_copy.applyMutations(change_list, chain_id)
    fixer_copy.findMissingResidues()
    fixer_copy.findMissingAtoms()
    fixer_copy.addMissingAtoms()
    fixer_copy.addMissingHydrogens()
    added_hydrogens = False
    error_count = 0
    success = True
    while not added_hydrogens:
#        print(fixer_test.positions[0][0]._value)
        if math.isnan(fixer_copy.positions[0][0]._value):
            print("failed to mutate residue...\nTrying again!\n")
            fixer_copy = copy.deepcopy(fixer)
            fixer_copy.applyMutations(change_list, chain_id)
            fixer_copy.findMissingResidues()
            fixer_copy.findMissingAtoms()
            fixer_copy.addMissingHydrogens()
            error_count += 1
            if error_count > max_iterations:
                success = False
                break
        else:
            added_hydrogens = True
    return fixer_copy, success

def apply_mutations(
        fixers, change_list, references, spring_k=15., bottom_width=0.005,
        rattle_dist=0.6, e_cutoff=-0.0):
    new_fixers = []
    energies = []
    min_attempts = 5
    success = True
    for num in range(len(fixers)):
        min_success = False
        attempts = 0
        while min_success != True:
            print("    mutating structure "+str(num))
            new_fixer,mut_success = apply_mutation(fixers[num],change_list)
            if mut_success:
                iis_struct = get_iis_to_constrain(
                    new_fixer, change_list, rattle_dist=rattle_dist)
                min_fixer,energy = get_fixer_energy(
                    new_fixer, references[num], spring_k=spring_k,
                    bottom_width=bottom_width, iis_struct=iis_struct)
                if (energy != None) and (np.isnan(energy) == False) and (energy < e_cutoff):
                    new_fixers.append(min_fixer)
                    energies.append(energy)
                    min_success = True
                else:
                    if attempts > min_attempts:
                        print("  Mutation FAILED!!")
                        break
                    else:
                        attempts += 1
                        print("  Retrying mutation! Attempt number "+str(attempts))
            else:
                break
        if (mut_success != True) or (min_success != True):
            success = False
            break
    return new_fixers,energies,success

def get_iis_to_constrain(fixer, change_list, rattle_dist=0.5):
    """Returns the indices that are not around the recent mutation to
       constrain with a flat-bottom potential. This function takes in
       a list of sequence changes but assumes that the list only contains
       a single mutation. Must change in the future to handle multiple
       mutations.
    """
    trj = fixer_to_trj(fixer)
    top = trj.topology
    res_number = int(change_list[0].split("-")[1])
    res_name = change_list[0].split("-")[2]
    # Get resid num and ensure mutation
    res_names = np.array(
        [str(res.name) for res in top.residues])
    # Residue numbering is wrong when converted directry
    # to mdtraj topology, so residue numbers are pulled from the fixer topology
    res_nums = np.array(
        [int(res.id) for res in fixer.topology.residues()])
    resid = np.where(res_nums==res_number)[0][0]
    if res_name != res_names[resid]:
        print(
            "SERIOUS ERROR!!!\nResidue "+str(res_number)+" in position "+\
            str(resid)+" does not match topology!")
    # Obtain indices that are not surrounding mutated residue
    res_iis = top.select("resi "+str(resid))
    neighbors = md.compute_neighbors(
        trj, cutoff=rattle_dist, query_indices=res_iis, periodic=False)[0]
    all_iis = top.select('all')
    constrain_iis = np.setdiff1d(all_iis,neighbors)
    return constrain_iis

def update_mdtraj_residues(fixer_top,mdtraj_top):
    res_names = np.array(
        [res.name for res in list(fixer_top.residues())],dtype=str)
    res_ids = np.array(
        [res.id for res in list(fixer_top.residues())],dtype=int)
    for num in range(len(list(mdtraj_top.residues))):
        list(mdtraj_top.residues)[num].name = res_names[num]
        list(mdtraj_top.residues)[num].resSeq = res_ids[num]
    return mdtraj_top

def save_fixers(fixers, filenames, output_directory=""):
    output_directory = os.path.abspath(output_directory)+"/"
    if not os.path.exists(output_directory):
        cmd = 'mkdir '+output_directory
        p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        out,err = p.communicate()
    for fixer_num in range(len(fixers)):
        try:
            fixer_top = fixers[fixer_num].topology
            mdtraj_top = md.Topology.from_openmm(fixer_top)
            updated_mdtraj_top = update_mdtraj_residues(fixer_top,mdtraj_top)
            fixer_pos = fixers[fixer_num].positions
            pdb = basics.create_mdtraj_from_pos(fixer_pos,updated_mdtraj_top)
            pdb.save_pdb(output_directory+filenames[fixer_num])
        except:
            print("write failed...")
    return

def fixer_to_trj(fixer):
    struct_top = md.Topology.from_openmm(fixer.topology)
    struct_pos = fixer.positions
    trj = basics.create_mdtraj_from_pos(struct_pos,struct_top)
    return trj

def relax_structure_old(
        fixer, pdb_ref, steps=1000, spring_k=15., bottom_width=0.005,
        iis_struct=None, iis_ref=None):
    """fixer is an openMM object, pdb_ref is an MDTraj object"""
    struct = fixer_to_trj(fixer)
    if iis_struct == None:
        iis_struct = struct.topology.select("backbone")
    if iis_ref == None:
        iis_ref = pdb_ref.topology.select("backbone")
    sim_flat_bottom = basics.setup_basic_sim_cartesian_restraint(
        struct, pdb_ref, spring_k, atm_inds=iis_struct,
        ref_atm_inds=iis_ref, bottom_width=bottom_width)
    sim_flat_bottom.context.setPositions(fixer.positions)
    sim_flat_bottom.step(int(steps))
    state = sim_flat_bottom.context.getState(getPositions=True)
    fixer.positions = state.getPositions(asNumpy=True)
    return fixer

def relax_structure(
        fixer, pdb_ref, steps=1000, spring_k=15., bottom_width=0.005,
        iis_struct=None):
    """fixer is an openMM object, pdb_ref is an MDTraj object"""
    struct = fixer_to_trj(fixer)
    sim_flat_bottom = basics.setup_sim_cartesian_restraint_backbone(
        struct, pdb_ref, spring_k, atm_inds=iis_struct,
        bottom_width=bottom_width)
    sim_flat_bottom.context.setPositions(fixer.positions)
    sim_flat_bottom.step(int(steps))
    state = sim_flat_bottom.context.getState(getPositions=True)
    fixer.positions = state.getPositions(asNumpy=True)
    return fixer

def get_fixer_energy(
        fixer,reference, max_iters1=200, max_iters2=2000, md_steps=1500,
        spring_k=15., bottom_width=0.005, iis_struct=None, iis_ref=None,
        verbose=False):
    t0 = time.time()
    if verbose:
        print("        -pre-minimizing structure")
    minimized = minimizers.minimize(fixer_to_trj(fixer),max_iterations=max_iters1,retry=False)
    min_fixer = copy.deepcopy(fixer)
    min_fixer.positions = minimized.openmm_positions(0)
    if verbose:
        print("        -relaxing structure")
    relaxed = relax_structure(
        min_fixer, reference, steps=md_steps, spring_k=spring_k,
        bottom_width=bottom_width, iis_struct=iis_struct)
    if verbose:
        print("        -post-minimizing structure")
    minimized = minimizers.minimize(fixer_to_trj(relaxed),max_iterations=max_iters2,retry=False)
    t1 = time.time()
    print("      openMM minimize took: "+str(t1-t0)+" s")
    energy = basics.get_energy(minimized)
    min_fixer = copy.deepcopy(fixer)
    min_fixer.positions = minimized.openmm_positions(0)
    return min_fixer,energy

def anneal_fixer_sidechains(
        fixer, reference, spring_k=15., bottom_width=0.005, max_T=350,
        min_T=200, T_spacing=0.1, steps_per_T=15):
    t0 = time.time()
    print("        -pre-minimizing structure")
    struct = fixer_to_trj(fixer)
    minimized = minimizers.minimize(struct, max_iterations=200,retry=False)
    # Setup anneal
    print("        -annealing structure")
    iis_struct = struct.topology.select("backbone")
    iis_ref = reference.topology.select("backbone")
    sim_flat_bottom = basics.setup_basic_sim_cartesian_restraint(
        minimized, reference, spring_k, atm_inds=iis_struct,
        ref_atm_inds=iis_ref, bottom_width=bottom_width)
    annealed = minimizers.simulated_annealing(
        minimized, sim=sim_flat_bottom, max_T=max_T, min_T=min_T,
        T_spacing=T_spacing, steps_per_T=steps_per_T)
    print("        -post-minimizing structure")
    minimized = minimizers.minimize(annealed,max_iterations=2000,retry=False)
    t1 = time.time()
    print("        openMM took: "+str(t1-t0)+" s")
    energy = basics.get_energy(minimized)
    final_fixer = copy.deepcopy(fixer)
    final_fixer.positions = minimized.openmm_positions(0)
    return final_fixer,energy 

def get_initial_fixer_energies(fixers,references):
    min_fixers = []
    energies = []
    for num in range(len(fixers)):
        print("    Getting initial energy for structure "+str(num))
        min_fixer,energy = anneal_fixer_sidechains(
            fixers[num], references[num], steps_per_T=50)
        min_fixers.append(min_fixer)
        energies.append(energy)
    return min_fixers,energies

def is_neighbor(seq1,seq2):
    if len(np.where(seq1!=seq2)[0])<2:
        neighbor = True
    else:
        neighbor = False
    return neighbor

def get_aij_from_neighbors(seqs):
    aij = np.zeros((len(seqs),len(seqs)),dtype=int)
    for i in range(len(seqs)):
        for j in range(i+1,len(seqs)):
            if is_neighbor(seqs[i],seqs[j]):
                aij[i,j] = 1
                aij[j,i] = 1
    return aij

def append_aij(aij, seqs):
    n_aij = len(aij)
    n_new_aij = len(seqs)
    if n_aij == n_new_aij:
        print("Aij matrix is up to date!")
    elif n_aij > n_new_aij:
        print("Error: Aij matrix has more entries than new sequences!")
        sys.exit()
    new_aij = np.zeros((n_new_aij,n_new_aij),dtype=int)
    new_aij[:n_aij,:n_aij] = aij
    for i in range(n_aij,n_new_aij):
        for j in range(n_aij):
            if is_neighbor(seqs[i],seqs[j]):
                new_aij[i,j] = 1
                new_aij[j,i] = 1
    for i in range(n_aij,n_new_aij):
        for j in range(i+1,n_new_aij):
            if is_neighbor(seqs[i],seqs[j]):
                new_aij[i,j] = 1
                new_aij[j,i] = 1
    return new_aij

def get_energy_ranks(energies):
    dev_sqs = []
    diffs = []
    for energy in energies:
        energy_slope,energy_intercept = np.polyfit(
            [0,len(energy)],[energy[0],energy[-1]],1)
        dev_sq = 0
        for num in range(1,len(energy)):
            dev_sq += ((energy_slope*num+energy_intercept)-energy[num])**2
        dev_sqs.append(dev_sq)
        diffs.append((energy[-1]-energy[0]))
    return np.array(diffs),np.array(dev_sqs)

def get_energy_ranks2(energies):
    """Assumes first half of energies are to be stabilized,
       and the second half to be destabilized
    """
    midpoint = int(len(energies[0])/2)
    diffs = []
    for energy in energies:
        diffs.append((-np.mean(energy[midpoint:])+np.mean(energy[:midpoint])))
    stabilities = np.mean(energies[:,:midpoint],axis=1)
    return np.array(diffs), np.array(stabilities)

def get_seq_diffs(seqs,ii=0):
    test_seq = seqs[ii]
    aa_diffs = []
    for seq in seqs:
        aa_diffs.append(len(np.where(test_seq!=seq)[0]))
    return np.array(aa_diffs)

def get_populations(energies,kT=100.):
    exps = np.exp(-energies/kT)
    pis = exps/(exps.sum(axis=1)[:,None])
    return pis

def average_energy(energies,kT):
    pis = get_populations(energies,kT)
    avg = np.sum(energies*pis,axis=1)
    return avg

def return_pops(energies, unfolded_energies, kT):
    state0 = average_energy(energies[:,:5],kT)
    state1 = average_energy(energies[:,5:10],kT)
    state2 = average_energy(energies[:,10:15],kT)
    state3 = average_energy(energies[:,15:20],kT)
    state4 = average_energy(energies[:,20:25],kT)
    state5 = average_energy(energies[:,25:],kT)
    unfolded_state = average_energy(unfolded_energies,kT)
    ensemble_list = np.array([unfolded_state,state0,state1,state2,state3,state4,state5]).T
    ensemble_pops = get_populations(ensemble_list,kT)
    unfolded_pops = ensemble_pops[:,0]
    closed_sites = ensemble_pops[:,1:4]
    open_sites = ensemble_pops[:,4:]
    return unfolded_pops,closed_sites,open_sites

def pick_test_seq_ii(seqs,energies,aij_connections,criterion='spreading',max_mut=4):
    if criterion == 'random':
        ii = np.random.choice(range(len(seqs)))
    elif criterion == 'spreading':
        aij = PageRank.Generate_Aij_Base(aij_connections)
        page_ranks = PageRank.Rank_Aij(aij)
        ii = np.argmin(page_rankings)
    elif criterion == 'FASTEnergy':
        aij = PageRank.Generate_Aij_Base(aij_connections)
        page_ranks = PageRank.Rank_Aij(aij,d=0.75)
        diff_ranks,stable_ranks = get_energy_ranks2(energies)
        diff_ranks_fs = FAST_Rankings.Feature_Scale(diff_ranks,minimize=True)
        stable_ranks_fs = FAST_Rankings.Feature_Scale(stable_ranks,minimize=True)
        page_ranks_fs = FAST_Rankings.Feature_Scale(page_ranks,minimize=True)
        seq_diffs = get_seq_diffs(seqs)
        seq_ranks = np.zeros(len(seqs))
        seq_ranks[np.where(seq_diffs<max_mut)] = 1
        alpha = 0.5
        beta = 0.25
        gamma = 0.25
        ii = np.argmax(
            (alpha*diff_ranks_fs)+
            (beta*stable_ranks_fs)+
            (gamma*page_ranks_fs)+
            seq_ranks)
    elif criterion == 'single':
        ii = 0
    return ii

def get_n_atoms(fixer):
    return len(list(fixer.topology.atoms()))

def format_energies(energies, digits=4):
    energy_list = ""
    for energy in energies:
        fmt_energy = ('{0:0.'+str(digits)+'f}').format(energy)
        energy_list += fmt_energy+" "
    return energy_list

def start_runs(
        output_directory, pdb_filenames, allowed_mutation_list,
        continue_run=False, unfolded_filenames=None):
    # Get directories and base things
    list_directory = output_directory+"Lists/"
    references = check_topologies(pdb_filenames)
    unfolded_references = check_topologies(unfolded_filenames)
    residue_list = get_residue_list(
        pdb_filenames[0], allowed_residues=allowed_mutation_list)
    pdb_filenames_states = [
        filename.split("/")[-1] for filename in pdb_filenames]
    unfolded_pdb_filenames_states = [
        filename.split("/")[-1] for filename in unfolded_filenames]
    if not continue_run:
        null_out = sp.check_output(["mkdir",output_directory])
        null_out = sp.check_output(["mkdir",list_directory])
        np.save(list_directory+"residue_list.npy",residue_list)
        seq_list = [copy.copy(residue_list['name'])]
        seq_list1d = ["".join(convert_3letter_code(residue_list['name']))]
        fixers = load_fixers(pdb_filenames)
        unfolded_fixers = load_fixers(unfolded_filenames)
        # Get initial energies
        print("\nAnnealing initial structures...")
        fixers,energies = get_initial_fixer_energies(fixers,references)
        print("\nAnnealing initial unfolded structures...")
        unfolded_fixers,unfolded_energies = get_initial_fixer_energies(
            unfolded_fixers,unfolded_references)
        n_atoms = [get_n_atoms(fixers[0])]
        energy_list = [energies]
        unfolded_energy_list = [unfolded_energies]
        Aij = [0]
        seq_ii = 0
        print("Initial energies = "+format_energies(energies))
        print("Initial unfolded energies = "+format_energies(unfolded_energies))
        print(
            "Initial energies/unfolded-avg = "+
            format_energies(energies/np.array(unfolded_energies).mean()))
        save_fixers(
            fixers, 0, pdb_filenames_states,
            output_directory=output_directory, new_directory=True)
        save_fixers(
            unfolded_fixers, 0, unfolded_pdb_filenames_states,
            output_directory=output_directory, new_directory=False)
        np.save(list_directory+"sequence_list_3letter.npy",seq_list)
        np.save(list_directory+"sequence_list.npy",seq_list1d)
        np.save(list_directory+"energy_list.npy",energy_list)
        np.save(list_directory+"unfolded_energy_list.npy",unfolded_energy_list)
        np.save(list_directory+"n_atoms_list.npy",n_atoms)
        np.save(list_directory+"Aij_Matrix.npy",Aij)
    else:
        print("\nContinuing from previous sequence search!\n")
        seq_list = np.load(list_directory+"sequence_list_3letter.npy")
        seq_list = [seq for seq in seq_list]
        seq_list1d = np.load(list_directory+"sequence_list.npy")
        seq_list1d = [seq for seq in seq_list1d]
        energy_list = np.load(list_directory+"energy_list.npy")
        energy_list = [energy for energy in energy_list]
        unfolded_energy_list = np.load(list_directory+"unfolded_energy_list.npy")
        unfolded_energy_list = [energy for energy in unfolded_energy_list]
        n_atoms = np.load(list_directory+"n_atoms_list.npy")
        n_atoms = [n_atom for n_atom in n_atoms]
        Aij = np.load(list_directory+"Aij_Matrix.npy")
        avg_unfolded = np.array(unfolded_energy_list).mean(axis=1)
        energies_per_unfolded = -np.array(energy_list)/avg_unfolded[:,None]
#        seq_ii = pick_test_seq_ii(
#            seq_list, energies_per_unfolded, Aij, criterion='FASTEnergy')
        seq_ii = 0
        fixers = load_fixers(
            pdb_filenames_states, output_directory=output_directory,
            mut_num=seq_ii)
        unfolded_fixers = load_fixers(
            unfolded_pdb_filenames_states, output_directory=output_directory,
            mut_num=seq_ii)
        residue_list['name'] = seq_list[seq_ii]
    start_num = len(seq_list)
    return \
        list_directory, references, unfolded_references, residue_list, \
        seq_list, seq_list1d, fixers, unfolded_fixers, n_atoms, energy_list, \
        unfolded_energy_list, Aij, seq_ii, start_num

###########################################################
# Main Function
###########################################################

def entry_point():

    t0 = time.time()
    print("\n##########   FELT   ##########")
    print("   Free Energy Landscape Tool\n")

    # Get residues to mutate
    allowed_mutation_list = np.loadtxt(
        "./Data/allowable_mutations.dat",dtype=int)

    # Get filenames from directory and load pdbs
    output_directory = "./Stabilize_State0_Attempt14/"
    pdb_paths = [
        "./Rep_States_Open_Closed/State0/",
        "./Rep_States_Open_Closed/State1228/",
        "./Rep_States_Open_Closed/State3315/",
        "./Rep_States_Open_Closed/State446/",
        "./Rep_States_Open_Closed/State614/",
        "./Rep_States_Open_Closed/State63/"]
    unfolded_pdb_paths = [
        "./RandomCoils/Extended_Structures/"]
    pdb_filenames = np.concatenate(
        [get_pdb_filenames(pdb_path) for pdb_path in pdb_paths])
    unfolded_pdb_filenames = np.concatenate(
        [get_pdb_filenames(pdb_path) for pdb_path in unfolded_pdb_paths])
    print("State-space filenames:")
    print(pdb_filenames)
    print("Unfolded filenames:")
    print(unfolded_pdb_filenames)
    pdb_filenames_states = [
        filename.split("/")[-1] for filename in pdb_filenames]
    unfolded_pdb_filenames_states = [
        filename.split("/")[-1] for filename in unfolded_pdb_filenames]

    continue_run = True
    list_directory, references, unfolded_references, residue_list, \
    seq_list, seq_list1d, fixers, unfolded_fixers, n_atoms, energy_list, \
    unfolded_energy_list, Aij, seq_ii, start_num = start_runs(
        output_directory, pdb_filenames, allowed_mutation_list,
        continue_run=continue_run, unfolded_filenames=unfolded_pdb_filenames)

    # Iterate mutations
    n_runs = 5000
    e_cutoff = -30000
    e_cutoff_unfolded = -20000
    for run_num in range(start_num,n_runs):
        success = False
        t_start = time.time()
        print("\nRUN"+str(run_num)+":")
        print("  Choosing mutation num "+str(seq_ii))
        attempts = 0
        while success != True:
            if attempts > 2:
                print("Too many failues! Something went wrong...")
                sys.exit()
            change_list, new_seq, new_seq_1d = generate_random_mutations(
                residue_list['number'], residue_list['name'], seq_list1d)
            t1 = time.time()
            new_fixers,energies,success1 = apply_mutations(
                fixers, change_list, references, rattle_dist=0.6, e_cutoff=e_cutoff)
            if success1:
                new_unfolded_fixers,unfolded_energies,success2 = apply_mutations(
                    unfolded_fixers, change_list, unfolded_references, spring_k=10.,
                    bottom_width=0.04, rattle_dist=0.6, e_cutoff=e_cutoff_unfolded)
                success = success2
            else:
                success = False
            if success:
                seq_list.append(new_seq)
                seq_list1d.append(new_seq_1d)
                Aij = append_aij(Aij,seq_list)
                t2 = time.time()
            attempts += 1
        n_atoms.append(get_n_atoms(new_fixers[0]))
        energy_list.append(energies)
        unfolded_energy_list.append(unfolded_energies)
        print("  Energies = "+format_energies(energies))
        print("  Unfolded-Energies = "+format_energies(unfolded_energies))
        print(
            "  Energies/Unfolded = "+
            format_energies(np.array(energies)/float(n_atoms[-1])))

        # Save lists and mutated PDBs
        save_fixers(
            new_fixers, run_num, pdb_filenames_states,
            output_directory=output_directory, new_directory=True)
        save_fixers(
            new_unfolded_fixers, run_num, unfolded_pdb_filenames_states,
            output_directory=output_directory, new_directory=False)
        np.save(list_directory+"sequence_list_3letter.npy",seq_list)
        np.save(list_directory+"sequence_list.npy",seq_list1d)
        np.save(list_directory+"energy_list.npy",energy_list)
        np.save(list_directory+"unfolded_energy_list.npy",unfolded_energy_list)
        np.save(list_directory+"n_atoms_list.npy",n_atoms)
        np.save(list_directory+"Aij_Matrix.npy",Aij)

        # Pick new state to mutate
        if run_num < 200:
            seq_ii = 0
        else:
            avg_unfolded = np.array(unfolded_energy_list).mean(axis=1)
            energies_per_unfolded = -np.array(energy_list)/avg_unfolded[:,None]
            seq_ii = pick_test_seq_ii(
                seq_list, energies_per_unfolded, Aij, criterion='FASTEnergy')
        residue_list['name'] = seq_list[seq_ii]
        fixers = load_fixers(
            pdb_filenames_states, output_directory=output_directory,
            mut_num=seq_ii)
        unfolded_fixers = load_fixers(
            unfolded_pdb_filenames_states, output_directory=output_directory,
            mut_num=seq_ii)
        t_end = time.time()
        print("  Run took "+str(t_end-t_start)+" s")
        print("    Generating ChangeList: "+str(t1-t_start)+" s")
#        print("    Appending Aij: "+str(t2-t1)+" s")
        print("    Applying Mutation: "+str(t2-t1)+" s")
        gc.collect()
    t4 = time.time()
    print("\nTotal time: "+str(t4-t0)+" s\n")

if __name__=='__main__':
    entry_point()
