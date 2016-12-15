# Author: Maxwell Zimmerman <mizimmer@wustl.edu>
# Copyright (c) 2016, Washington University
# All rights reserved.

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from __future__ import absolute_import, print_function, division
import numpy as np
import os
import sys
import time
from ..input_control import loading

__all__=['core_sampling']

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

class core_sampling():
    def __init__(self):
        self.structure_filenames = self.args.structure_filenames
        self.runs = self.args.runs
        self.max_mutations = self.args.max_mutations
        self.residues_and_mutations = self.args.residues_and_mutations
        self.output_name = self.args.output_name
        self.forcefield = self.args.forcefield
        self.energy_error = self.args.energy_error
        self.spring_const = self.args.spring_const
        self.bottom_width = self.args.bottom_width
        self.rattle_distance = self.args.rattle_distance
        self.simulation_steps = self.args.simulation_steps
        self.postmin_steps = self.args.postmin_steps
        self.anneal = self.args.anneal
        self.anneal_spring_const = self.args.anneal_spring_const
        self.anneal_temp_range = self.anneal_temp_range
        self.anneal_steps = self.args.anneal_steps
        self.rotate_chi1 = self.args.rotate_chi1

    def print_basic_inputs(self):
        print("Structure filenames:\t\t"+str(self.structure_filenames))
        print("Number of runs:\t\t\t"+str(self.runs))
        print("Maximum sequence mutations:\t"+str(self.max_mutations))
        print("Allowable mutations:\t\t"+str(self.residues_and_mutations))
        print("Output directory:\t\t"+str(self.output_name))
        print("Protein Forcefield:\t\t"+str(self.forcefield))
        print("Energy error cutoff:\t\t"+str(self.energy_error)+" kJ/mol")
        print("MD spring const:\t\t"+str(self.spring_const)+" kJ/mol/nm")
        print("MD bottomw width:\t\t"+str(self.bottom_width)+" nm")
        print("MD rattle distance:\t\t"+str(self.rattle_distance)+" nm")
        print("MD steps:\t\t\t"+str(self.simulation_steps))
        print("Post-minimization steps:\t"+str(self.postmin_steps))
        if self.anneal:
            print("Annealing initial sidechains:\tTrue")
            print("Annealing spring_const:\t\t"+str(self.anneal_spring_const)+" kJ/mol/nm")
            print("Annealing temperature range:\t"+str(self.anneal_temp_range)+" K")
            print("Anneal MD steps:\t\t"+str(self.anneal_steps))
        else:
            print("Annealing initial sidechains:\tFalse")
        print("Rotate mutations by chi1:\t"+str(self.rotate_chi1))

    def check_basic_inputs(self):
        self.structure_filenames = loading.check_filenames(
            self.structure_filenames)
        loading.check_runs(self.runs)
        loading.check_max_mutations(self.max_mutations)
        self.residues_and_mutations = \
            loading.check_residues_and_mutations(
                self.residues_and_mutatons)









