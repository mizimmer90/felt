# Author: Maxwell Zimmerman <mizimmer@wustl.edu>
# Copyright (c) 2016, Washington University
# All rights reserved.

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import numpy as np
import os
import sys
from ..exceptions import ImproperStructureFiles, ImproperInputRange, ImproperMutations
from ..mutation_tools import pdb_tools

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------


def check_filenames(filenames):
    if len(filenames) == 1:
        try:
            struct_file = open(filenames[0],'r')
            filenames = [
                filename.split()[0] for filename in struct_file.readlines()]
            struct_file.close()
        except:
            raise ImproperStructureFiles(
                "Only one file detected, and it is"+\
                "not a text file of more filenames.")
    for filename in filenames:
        if type(filename) is not str:
            raise ImproperStructureFiles(
                "Type {} is not string.".format(type(filename)))
        if filename.split(".")[-1] != 'pdb':
            raise ImproperStructureFiles(
                "The file '{}' extension is not a pdb.".format(filename))
        if not os.path.exists(filename):
            raise ImproperStructureFiles(
                "The file '{}' does not exist.".format(filename))
    return filenames

def check_runs(runs):
    if runs <= 0:
        raise ImproperInputRange(
            'Number of runs must be greater than zero.')
    return

def check_max_mutations(max_mutations):
    if max_mutations <= 0:
        raise ImproperInputRange(
            'Maximum number of mutations must be greater than zero.')
    return

def check_residues_and_mutations(residues_and_mutations):
    if not os.path.exists(residus_and_mutations):
        raise ImproperMutations(
            'The file "{}" does not exist'.format(residues_and_mutations))
    res_file = open(residues_and_mutations,'r')
    res_file_info = [line.split() for line in res_file.readlines()]
    res_file.close()
    res_data = []
    for line_num in range(len(res_file_info)):
        res_num = int(res_file_info[0])
        if len(res_file_info[line_num]) > 3:
            raise ImproperMutations(
                'Unrecognized mutation format. Should be '+\
                '{res_number (not) (allowed_residues)}\n'+\
                'i.e. "10 not KLPT" or "10 AGHRTY"')
        elif len(res_file_info[line_num]) == 2:
            allowed_aas = res_file_info[line_num][1]

    return














