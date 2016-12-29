# Author: Maxwell I. Zimmerman
# Copyright (c) 2016, Washington University
# Washington University in St. Louis

##############################################################################
# Imports
##############################################################################


from __future__ import absolute_import, print_function, division
import copy
import numpy as np
import os
import sys
import time
from ..exceptions import PathExists,UnexpectedError
from ..input_output import loading, output
from ..tools import minimizers,pdb_tools,sim_basics,utils

__all__=['core_sampling']

##############################################################################
# Code
##############################################################################

class core_sampling(object):
    def __init__(self):
        # Initialize arguments
        self.structure_filenames = self.args.structure_filenames
        self.runs = self.args.runs
        self.max_mutations = self.args.max_mutations
        self.residues_and_mutations = self.args.residues_and_mutations
        self.output_directory = self.args.output_directory
        self.continue_sampling = self.args.continue_sampling
        self.forcefield = self.args.forcefield
        self.sol_forcefield = self.args.sol_forcefield
        self.energy_error = self.args.energy_error
        self.spring_const = self.args.spring_const
        self.bottom_width = self.args.bottom_width
        self.rattle_distance = self.args.rattle_distance
        self.simulation_steps = self.args.simulation_steps
        self.postmin_steps = self.args.postmin_steps
        self.anneal = self.args.anneal
        self.anneal_spring_const = self.args.anneal_spring_const
        self.anneal_temp_range = self.args.anneal_temp_range
        self.anneal_steps = self.args.anneal_steps
        self.rotate_chi1 = self.args.rotate_chi1

    def print_basic_inputs(self):
        print("Structure filenames:\t\t"+str(self.structure_filenames))
        print("Number of runs:\t\t\t"+str(self.runs))
        print("Maximum sequence mutations:\t"+str(self.max_mutations))
        print("Allowable mutations:\t\t"+str(self.residues_and_mutations))
        if self.output_directory is None:
            self.output_directory = self.klass.__name__
        print("Output directory:\t\t"+str(self.output_directory))
        print("Protein Forcefield:\t\t"+str(self.forcefield))
        print("Solvation Forcefield:\t\t"+str(self.sol_forcefield))
        print("Energy error cutoff:\t\t"+str(self.energy_error)+" kJ/mol")
        print("MD spring const:\t\t"+str(self.spring_const)+" kJ/mol/nm")
        print("MD bottom width:\t\t"+str(self.bottom_width)+" nm")
        print("MD rattle distance:\t\t"+str(self.rattle_distance)+" nm")
        print("MD steps:\t\t\t"+str(self.simulation_steps))
        print("Post-minimization steps:\t"+str(self.postmin_steps))
        if self.anneal == 'True':
            print("Annealing initial sidechains:\tTrue")
            print("Annealing spring_const:\t\t"+str(self.anneal_spring_const)+" kJ/mol/nm")
            print("Annealing temperature range:\t"+str(self.anneal_temp_range)+" K")
            print("Anneal MD steps:\t\t"+str(self.anneal_steps))
        else:
            print("Annealing initial sidechains:\tFalse")
        print("Rotate mutations by chi1:\t"+str(self.rotate_chi1))

    def load_basic_inputs(self):
        # Print time and program status update
        output.output_status('loading input data')
        # Check the validity of structure filenames and load mdtraj 
        # and PDBFixer objects.
        self.structure_filenames = loading.check_filenames(
            self.structure_filenames)
        self.base_filenames = [
            filename.split("/")[-1] for filename in self.structure_filenames]
        self.references = loading.load_references(self.structure_filenames)
        self.fixers = loading.load_fixers(self.structure_filenames)
        # Check that the runs specified and maximum mutations are reasonable.
        loading.check_runs(self.runs)
        loading.check_max_mutations(self.max_mutations)
        # If an allowable mutation list is specified, load it and ensure that
        # amino acids are valid.
        if self.residues_and_mutations:
            self.residues_and_mutations = \
                loading.check_residues_and_mutations(
                    self.residues_and_mutations)
        else:
            self.residues_and_mutations = \
                loading.generate_residues_and_mutations(
                    self.references[0])
        # Check the range of various MD parameters
        loading.check_spring_const(self.spring_const)
        loading.check_bottom_width(self.bottom_width)
        loading.check_rattle_distance(self.rattle_distance)
        loading.check_simulation_steps(self.simulation_steps)
        loading.check_postmin_steps(self.postmin_steps)
        if self.anneal == 'True':
            loading.check_anneal_spring_const(self.anneal_spring_const)
            loading.check_anneal_steps(self.anneal_steps)
            self.anneal_temp_range = loading.check_anneal_temp_range(
                self.anneal_temp_range)
        # Update naming on forcefields
        self.forcefield += '.xml'
        self.sol_forcefield += '.xml'

    def create_directory_structure(self):
        # Get full path for output directory
        self.output_directory = os.path.abspath(self.output_directory)+"/"
        # Print time and program status update
        output.output_status(
            'creating directory structure: "%s"' % self.output_directory)
        # Check that main directory does not exist and create it
        if os.path.exists(self.output_directory):
            raise PathExists(
                'The specified output directory already exists! '+\
                'If restarting a previous run, specify the '+\
                'continue_sampling flag.')
        self.data_directory = self.output_directory+'Data/'
        cmd1 = 'mkdir '+self.output_directory
        cmd2 = 'mkdir '+self.data_directory
        cmds = [cmd1, cmd2]
        utils.run_commands(cmds)

    def initialize_variables(self):
        self.Aij = np.array([[0]])
        self.Aij_filename = os.path.abspath(self.data_directory+"Aij.npy")
        self.energies = []
        self.energies_filename = os.path.abspath(
            self.data_directory+"energies.npy")
        self.seqs_1d = []
        self.seqs_1d_filename = os.path.abspath(
            self.data_directory+"seqs_1d.npy")
        self.seqs_3letter = []
        self.seqs_3letter_filename = os.path.abspath(
            self.data_directory+"seqs_3letter.npy")
        self.run_to_mutate = 0

    def update_sampling_data(self):
        # Get full path for output directory
        self.output_directory = os.path.abspath(self.output_directory)+"/"
        self.data_directory = self.output_directory+'Data/'
        # Initialize and load important variables
        self.Aij_filename = os.path.abspath(self.data_directory+"Aij.npy")
        self.Aij = np.load(self.Aij_filename)
        self.energies_filename = os.path.abspath(
            self.data_directory+"energies.npy")
        self.energies = [
            list(energy) for energy in np.load(self.energies_filename)]
        self.seqs_1d_filename = os.path.abspath(
            self.data_directory+"seqs_1d.npy")
        self.seqs_1d = [
            seq for seq in np.load(self.seqs_1d_filename)]
        self.seqs_3letter_filename = os.path.abspath(
            self.data_directory+"seqs_3letter.npy")
        self.seqs_3letter = [
            list(seq) for seq in np.load(self.seqs_3letter_filename)]
        self.run_number = len(self.seqs_1d)-1
        self.run_directory = self.output_directory+'RUN'+\
            str(self.run_number)+'/'
        loading.check_important_variables(
            self.Aij,self.energies,self.seqs_1d,self.seqs_3letter)
        loading.check_run_numbers(self.output_directory,self.run_number+1)

    def anneal_initial_structures(self):
        output.output_status('annealing %d starting structures' % len(self.fixers)) 
        self.energies = [[0 for i in range(len(self.fixers))]]
        for fixer_num in range(len(self.fixers)):
            output.output_status('annealing structure %d' % fixer_num)
            self.fixers[fixer_num],self.energies[0][fixer_num] = \
                minimizers.anneal_fixer_sidechains(
                    self.fixers[fixer_num], spring_const=self.anneal_spring_const,
                    bottom_width=0.001, T_min=self.anneal_temp_range[0],
                    T_spacing=self.anneal_temp_range[1],
                    T_max=self.anneal_temp_range[2], steps_per_T=self.anneal_steps,
                    prot_ff=self.forcefield, sol_ff=self.sol_forcefield)

    def append_energies(self):
        new_energies = []
        for num in range(len(self.fixers)):
            pdb = sim_basics.pdb_from_fixer(self.fixers[num])
            sim = sim_basics.setup_basic_sim(
                pdb, prot_ff=self.forcefield, sol_ff=self.sol_forcefield)
            new_energies.append(sim_basics.get_energy(pdb, sim=sim))
        self.energies.append(new_energies)

    def append_seqs(self):
        new_seq_3letter = pdb_tools.get_sequence(
            self.fixers[0], res_subset=self.residues_and_mutations['res'])
        self.seqs_3letter.append(new_seq_3letter)
        self.seqs_1d.append(
            pdb_tools.convert_3letter_seq(
                new_seq_3letter, concat_output=True))

    def save_base_run_data(self):
        self.run_directory = self.output_directory+'RUN'+\
            str(self.run_number)+'/'
        cmd = 'mkdir %s' % self.run_directory
        utils.run_commands([cmd])
        # Save structures
        output.output_status('saving pdbs')
        output_filenames = [
            os.path.abspath(
                self.run_directory+filename)
            for filename in self.base_filenames]
        output.save_fixers_as_pdbs(self.fixers, output_names=output_filenames)
        # Get and save energies
#        output.output_status('saving energies')
#        core_sampling.append_energies(self)
        np.save(self.energies_filename,self.energies)
        # Get and save sequence info
        output.output_status('updating sequences')
        core_sampling.append_seqs(self)
        np.save(self.seqs_1d_filename, self.seqs_1d)
        np.save(self.seqs_3letter_filename, self.seqs_3letter)
        # Save adjacency matrix
        if self.run_number != 0:
            self.Aij = pdb_tools.append_aij(self.Aij, self.seqs_3letter)
        np.save(self.Aij_filename,self.Aij)

    def load_new_fixers(self, run_to_mutate):
        output.output_status('loading new pdbs to mutate')
        self.run_to_mutate = run_to_mutate
        self.run_directory = self.output_directory+'RUN'+\
            str(self.run_to_mutate)+'/'
        input_filenames = [
            os.path.abspath(
                self.run_directory+filename)
            for filename in self.base_filenames]
        self.fixers = loading.load_fixers(input_filenames)

    def select_new_mutations(self):
        unique_mutation_attempt = 0
        seq = self.seqs_3letter[self.run_to_mutate]
        while True:
            res_num,allowed_muts = np.random.choice(self.residues_and_mutations)
            res_ii = np.where(self.residues_and_mutations['res']==res_num)
            prev_res = self.seqs_3letter[self.run_to_mutate][res_ii[0][0]]
            mutation_list = pdb_tools.convert_1letter_seq(allowed_muts)
            new_res = pdb_tools.select_random_mutation(
                res_list=mutation_list, exclude=prev_res)
            # Test if new sequence is unique
            new_seq = copy.copy(seq)
            new_seq[res_ii[0][0]] = new_res
            new_seq_1d = pdb_tools.convert_3letter_seq(
                new_seq,concat_output=True)
            if not np.any(np.array(self.seqs_1d)==new_seq_1d):
                break
            unique_mutation_attempt += 1
            if unique_mutation_attempt >= 200:
                raise UnexpectedError(
                    'Unable to discover a mutation that generates a '+\
                    'unique sequence!')
        output.output_mutation(res_num,prev_res,new_res)
        self.change_list = ['%s-%d-%s' % (prev_res, res_num, new_res)]

    def apply_mutations_to_fixers(self):
        new_fixers = []
        new_energies = []
        for num in range(len(self.fixers)):
            output.output_status('mutating structure %d' % num)
            mutated_fixer,success = pdb_tools._apply_mutations(
                self.fixers[num], self.change_list)
            pdb = sim_basics.pdb_from_fixer(mutated_fixer)
            restrain_iis = pdb_tools._get_restraint_iis(
                self.change_list, pdb=pdb,
                rattle_distance=self.rattle_distance)
            if  self.rotate_chi1:
                res_num = int(self.change_list[0].split("-")[1])
                res_name = self.change_list[0][-3:]
                if res_name != 'ALA' and res_name != 'PRO' \
                        and res_name != 'GLY':
                    output.output_status('rotating chi1 for better sampling')
                    pdbs = pdb_tools.rotate_chi1(pdb,res_num)
                else:
                    pdbs = [pdb]
            else:
                pdbs = [pdb]
            tmp_energies = []
            tmp_fixers = []
            for num in range(len(pdbs)):
                mutated_fixer.positions = pdbs[num].openmm_positions(0)
                rattled_fixer,energy = minimizers.relax_mutation(
                    mutated_fixer, self.references[num], max_iters1=0,
                    max_iters2=self.postmin_steps, md_steps=self.simulation_steps,
                    spring_const=self.spring_const, bottom_width=self.bottom_width,
                    iis_struct=restrain_iis, prot_ff=self.forcefield,
                    sol_ff=self.sol_forcefield)
                tmp_energies.append(energy)
                tmp_fixers.append(rattled_fixer)
            ii = np.argmin(tmp_energies)
            new_fixers.append(tmp_fixers[ii])
            new_energies.append(tmp_energies[ii])
        self.energies.append(new_energies)
        self.fixers = new_fixers

