# Author: Maxwell I. Zimmerman
# Copyright (c) 2016, Washington University
# Washington University in St. Louis

##############################################################################
# Imports
##############################################################################

import numpy as np
import time
from ..cmdline import NumpydocClassCommand, Command, argument, argument_group, exttype
from ..input_output import output
from ..sampling.FASTMutate import FASTMutate
from ..sampling.RandomMutate import RandomMutate

##############################################################################
# Imports
##############################################################################

class MutationSearch(NumpydocClassCommand):
    _group = '1-search'

    # Mutational stuff
    g1 = argument_group('Mutational Input')
    structure_filenames = g1.add_argument(
        '-i', '--structure_filenames', required=True, nargs='+',
        help='Path to the structure filenames or text file containing all filenames')
    runs = g1.add_argument(
        '-r', '--runs', required=False, default=100,
        type=int, help='Number of total mutations to attempt')
    max_mutations = g1.add_argument(
        '-m', '--max_mutations', required=False, type=int, default=np.inf,
        help='The run number to continue sampling from')
    residues_and_mutations = g1.add_argument(
        '-M', '--residues_and_mutations', required=False, type=str, default=None,
        help='File containing residue numbers and allowable mutations')
    output_directory = g1.add_argument(
        '-o', '--output_directory', required=False, default=None, type=str,
        help='The name for the directory output')
    continue_sampling =  g1.add_argument(
        '-c','--continue_sampling', required=False, default=False, action='store_true',
        help='Flag to continue sampling from previous run.')

    # Modeling stuff
    g2 = argument_group('Modeling Input')
    forcefield = g2.add_argument(
        '-ff', '--forcefield', required=False, choices=[
            'amber96','amber99sb','amber03','amber10',
            'amoeba2013','charmm_polar_2013'],
        default='amber03', help='''The forcefield to use for simulation and
            energy calculations''')
    sol_forcefield = g2.add_argument(
        '-sf', '--sol_forcefield', required=False, choices=['amber03_obc'],
        default='amber03_obc', help='''The solvent forcefield to use for
            simulation and energy calculations''')
    energy_error = g2.add_argument(
        '-ee', '--energy_error', required=False, default=0.0,
        help='''The maximum energy to be considered for a successful mutation.
            i.e. if an energy is calculated above this value, FELT will not
            consider the mutation in calculations. This is useful because
            errored mutations can skew population calculations. (kcal/mol)''')
    spring_const = g2.add_argument(
        '-k', '--spring_const', default=15.0, type=float, required=False,
        help='''The spring constant to use for restraining backbone and
            sidechains that are spatially proximate to the mutation site
            with a flat-welled potential.''')
    bottom_width = g2.add_argument(
        '-bw', '--bottom_width', required=False, type=float, default=0.005,
        help='''The bottom width of the flat-welled potential that is used
            on atoms spatially proximate to the mutation site''')
    rattle_distance = g2.add_argument(
        '-d', '--rattle_distance', required=False, type=float, default=0.1,
        help='''The distance from any atom in the mutated residue to apply a
            flat-welled potential to. Other residues are fixed.''')
    simulation_steps = g2.add_argument(
        '-s','--simulation_steps', required=False, type=int, default=1000,
        help='''The number of steps to perform of MD simulation to rattle
            sidechains.''')
    postmin_steps = g2.add_argument(
        '-P','--postmin_steps', required=False, type=int, default=1000,
        help='The number of minimization steps to perform after simulation.')
    anneal = g2.add_argument(
        '-a', '--anneal', required=False, choices=['True','False'],
        default='True', help='Optionally anneal sidechains of starting states.')
    anneal_spring_const = g2.add_argument(
        '-ak','--anneal_spring_const', required=False, type=float, default=15.0,
        help='''The spring constant to use for restraining backbone atoms 
            during the initial annealing steps.''')
    anneal_temp_range = g2.add_argument(
        '-at','--anneal_temp_range', required=False, type=str,
        default="200:0.1:350", help='''The temperature range to anneal
            sidechains if annealing initial structures. An input of 
            200:0.1:350 will start at 200K and increase the temp to
            350K in increments of 0.1K.''')
    anneal_steps = g2.add_argument(
        '-as', '--anneal_steps', required=False, type=int, default=15,
        help='''The number of dynamics steps to take at each temperature
            interval when annealing.''')
    rotate_chi1 = g2.add_argument(
        '-R', '--rotate_chi1', required=False, type=bool, default=False,
        help='''If selected, will rotate mutations around chi1 to get broader
            sampling. This generates 3x more structures to calculate energies.''')

    def __init__(self, args):
        self.args = args
        self.klass.__init__(self)
        print("\n                        Running " + self.klass.__name__)
        print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("~~~~                          Inputs                            ~~~~")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")


    def start(self):
        self.klass.print_and_check_inputs(self)
        if not self.continue_sampling:
            self.klass.create_directory_structure(self)
            if self.anneal:
                self.klass.anneal_initial_structures(self)
        else:
            output.output_status('continuing from previous run')
            self.klass.update_sampling_data(self)

class FASTMutateCommand(MutationSearch):
    klass = FASTMutate
    _concrete = True
    _group = '1-search'

class RandomMutateCommand(MutationSearch):
    klass = RandomMutate
    _concrete = True
    _group = '1-search'

