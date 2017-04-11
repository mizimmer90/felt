# Author: Maxwell I. Zimmerman
# Copyright (c) 2016, Washington University
# Washington University in St. Louis

##############################################################################
# Imports
##############################################################################

import glob
import mdtraj as md
import numpy as np
import os
import sys
from ..exceptions import ImproperStructureFiles, ImproperInputRange, ImproperMutations, InvalidData
from ..tools import pdb_tools
from pdbfixer import PDBFixer

##############################################################################
# Code
##############################################################################

def _load_data(filename):
    file_to_load = open(filename,'r')
    file_data = np.array(file_to_load.readlines(),dtype=str)
    file_to_load.close()
    return file_data[np.where(file_data!='\n')]

def check_filenames(filenames):
    if len(filenames) == 1:
        try:
            file_data = _load_data(filenames[0])
            filenames = [
                filename.split()[0] for filename in file_data]
        except:
            raise ImproperStructureFiles(
                "Only one file detected, and it is "+\
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

def load_references(filenames):
    try:
        pdbs = md.load(filenames[0])
    except:
        raise ImproperStructureFiles(
            'Could not load pdb file "{}"'.format(filenames[0]))
    for filename in filenames:
        try:
            pdbs = pdbs.join(md.load(filename))
        except:
            raise ImproperStructureFiles(
                'Could not load pdb file "{}"'.format(filename))
    return pdbs

def load_fixers(filenames):
    fixers = []
    for filename in filenames:
        fixers.append(PDBFixer(filename=filename))
    return fixers

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
    res_file_info = [
        line.split() for line in _load_data(residues_and_mutations)]
    res_data = []
    for line_num in range(len(res_file_info)):
        res_num = int(res_file_info[line_num][0])
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
        res_data.append((res_num,allowed_aas))
    res_data = np.array(res_data,dtype=[('res','int'),('seq',np.str_,20)])
    return res_data

def generate_residues_and_mutations(pdb):
    '''This generates an array of residues and allowed mutations from a pdb.
       It specifies that all residues can use all amino acids.'''
    allowed_aas = "".join(pdb_tools.protein_residues_1letter)
    res_data = np.array(
        [(res.resSeq,allowed_aas) for res in pdb.topology.residues],
        dtype=[('res','int'),('seq',np.str_,20)])
    return res_data

def check_spring_const(spring_const):
    if spring_const <= 0:
        raise ImproperInputRange(
            'Spring constant must be greater than zero.')
    return

def check_bottom_width(bottom_width):
    if bottom_width <= 0:
        raise ImproperInputRange(
            'Bottom width for flat-bottom potential must'+\
            ' be greater than zero.')
    return

def check_rattle_distance(rattle_distance):
    if rattle_distance <= 0:
        raise ImproperInputRange(
            'Mutation rattle distance must be greater than zero.')
    return

def check_simulation_steps(simulation_steps):
    if simulation_steps <= 0:
        raise ImproperInputRange(
            'Number of simulation steps must be greater than zero.')
    return

def check_postmin_steps(postmin_steps):
    if postmin_steps <= 0:
        raise ImproperInputRange(
            'The post-simulation minimization steps must be greater'+\
            ' than zero.')
    return

def check_anneal_spring_const(anneal_spring_const):
    if anneal_spring_const  <= 0:
        raise ImproperInputRange(
            'The spring constant for annealing must be greater than zero.')
    return

def check_anneal_steps(anneal_steps):
    if anneal_steps <= 0:
        raise ImproperInputRange(
            'The number of simulation steps per temperature interval when'+\
            ' annealing must be greater than zero.')
    return

def check_anneal_temp_range(temp_range):
    try:
        start,step,stop = [float(t) for t in temp_range.split(":")]
    except:
        raise ImproperInputRange(
            'The annealing temperature range format was not valid. '+\
            'A valid format would be "start:step:stop"')
    if start <= 0:
        raise ImproperInputRange(
            'The starting temperature for '+\
            'annealing must be greater than zero.')
    if step <= 0:
        raise ImproperInputRange(
            'The annealing step must be greater than zero.')
    if stop <= 0:
        raise ImproperInputRange(
            'The final temperature for '+\
            ' annealing must be greater than zero.')
    if stop <= start:
        raise ImproperInputRange(
            'The final annealing temperature must be greater '+\
            'than the starting temperature.')
    if step > (stop-start):
        raise ImproperInputRange(
            'The annealing step should be smaller than the difference '+\
            'between final and initial temperatures.')
    return [start,step,stop]

def check_run_numbers(directory,correct_nums):
    directories = glob.glob(directory+"RUN*")
    nums = np.array([int(path.split("RUN")[-1]) for path in directories])
    if len(nums) != correct_nums:
        raise ImproperInputRange(
            'The number of runs found (%d) does not match expected number of runs (%d)' \
            % (len(nums), correct_nums))
    if not np.all(nums==np.array(range(len(nums)))):
        raise ImproperInputRange('The provided runs are not sequential')
    return

def check_important_variables(Aij,energies,seqs_1d,seqs_3letter):
    Aij_first,Aij_second = Aij.shape
    if Aij_first != Aij_second:
        raise InvalidData('Aij matrix dimensions do not match.')
    energy_len = len(energies)
    seqs_1d_len = len(seqs_1d)
    seqs_3letter_len = len(seqs_3letter)
    if not ((Aij_first == energy_len) and (Aij_first == seqs_1d_len) and \
            (Aij_first == seqs_3letter_len)):
        raise InvalidData(
            'The number of entries in core variables are not consistent!\n'+\
            'Aij dimensions: %d\n' % Aij_first +\
            'energies: %d\n' % energy_len +\
            'sequence list 1d: %d\n' % seqs_1d_len+\
            'sequence list 3-letter: %d' % seqs_3letter_len)
    return

