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
from ..tools import pdb_tools

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

def _negate_mutations(not_allowed_aas):
    not_allowed_aas = [aa for aa in not_allowed_aas]
    allowed_aas = "".join(
        np.setdiff1d(pdb_tools.protein_residues_1letter, not_allowed_aas))
    return allowed_aas

def _test_aas(aas):
    recognized_aas = np.array(pdb_tools.protein_residues_1letter)
    for aa in aas:
        if not np.any(recognized_aas==aa):
            raise ImproperMutations(
                'Amino acid letter "{}" is not recognized'.format(aa))
    return

def check_residues_and_mutations(residues_and_mutations):
    if not os.path.exists(residues_and_mutations):
        raise ImproperMutations(
            'The file "{}" does not exist'.format(residues_and_mutations))
    res_file = open(residues_and_mutations,'r')
    res_file_info = [line.split() for line in res_file.readlines()]
    res_file.close()
    res_data = []
    for line_num in range(len(res_file_info)):
        res_num = int(res_file_info[line_num][0])
        print(res_num)
        if len(res_file_info[line_num]) == 1:
            allowed_aas = "".join(pdb_tools.protein_residues_1letter)
        elif len(res_file_info[line_num]) == 2:
            allowed_aas = res_file_info[line_num][1]
            _test_aas(allowed_aas)
        elif len(res_file_info[line_num]) == 3:
            negation = res_file_info[line_num][1]
            if negation != 'not':
                raise ImproperMutations(
                    'Second input on mutations line should either be "not"'+\
                    ' or a list of amino acids.')
            not_allowed_aas = res_file_info[line_num][2]
            _test_aas(not_allowed_aas)
            allowed_aas = _negate_mutations(not_allowed_aas)
        elif len(res_file_info[line_num]) > 3:
            raise ImproperMutations(
                'Unrecognized mutation format. Should be '+\
                '{res_number (not) (allowed_residues)}\n'+\
                'i.e. at position 10: \n'+\
                '    "10" for all residues\n'+
                '    "10 not KLPT" for all residues except KLPT\n'+
                '    "10 AGHRTY" for only residues AGHRTY\n')
        res_data.append([res_num,allowed_aas])
    return res_data












