# Author: Maxwell I. Zimmerman
# Copyright (c) 2016, Washington University
# Washington University in St. Louis

##############################################################################
# Imports
##############################################################################


from __future__ import absolute_import, print_function, division
import numpy as np
import os
import sys
import time
from ..exceptions import PathExists
from ..input_output import loading, output
from ..tools import minimizers,utils

__all__=['core_sampling']

##############################################################################
# Code
##############################################################################

class core_sampling(object):
    def __init__(self):
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
        else:
            self.run_number = 0
            self.run_directory = self.output_directory+'RUN'+\
                str(self.run_number)+'/'
            self.data_directory = self.output_directory+'Data/'
            cmd1 = 'mkdir '+self.output_directory
            cmd2 = 'mkdir '+self.run_directory
            cmd3 = 'mkdir '+self.data_directory
            cmds = [cmd1, cmd2, cmd3]
            utils.run_commands(cmds)

    def anneal_initial_structures(self):
        output.output_status('annealing starting structures')        
        for fixer_num in range(len(self.fixers)):
            self.fixers[fixer_num] = minimizers.anneal_fixer_sidechains(
                self.fixers[fixer_num], spring_const=self.anneal_spring_const,
                bottom_width=0.001, T_min=self.anneal_temp_range[0],
                T_spacing=self.anneal_temp_range[1],
                T_max=self.anneal_temp_range[2], steps_per_T=self.anneal_steps,
                prot_ff=self.forcefield, sol_ff=self.sol_forcefield)


    def update_sampling_data(self):
        pass










