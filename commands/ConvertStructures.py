"""Converts multiple gromacs files (pdb, gro, trr, xtc etc.).
"""
# Author: Maxwell Zimmerman <mizimmer@wustl.edu>
# Contributors:
# Copyright (c) 2016, Washington University
# All rights reserved.

import mdtraj as md
import numpy as np
import os
import time
from ..cmdline import Command, argument, argument_group, exttype
from ..utilities import tools

class ConvertStructures(Command):
    _group = '0-Support'
    _concrete = True
    description = __doc__
    pdbs = argument('-p', '--pdbs', nargs='+', required=False, default=None,
        help='Path to PDB files')
    gros = argument('-g', '--gros', nargs='+', required=False, default=None,
        help='Path to GRO files')
    structure = argument('-s','--structure', type=str, required=False,
        default='md.tpr', help='''
            The structure file to do conversions (required for trjconv)''')
    trajectories = argument('-t','--trajectories', nargs='+', required=False,
        default=None, help='List of trajectories to process')
    trjconv = argument('-m','--trjconv', choices=['backbone','masses','system'],
        default='system', help='''
            Does trjconv on multiple structures in parallel to center, remove
            pbc, and compact protein. Outputs either backbone, protein masses,
            or the system atomic coordinates.''')
    cores = argument('-c', '--cores', type=int, required=False, default=1,
        help='The number of processors to use for conversions')

    def __init__(self, args):
        self.args = args

    def start(self):
        pdbs = self.args.pdbs
        gros = self.args.gros
        structure = self.args.structure
        trajectories = self.args.trajectories
        trjconv = self.args.trjconv
        cores = self.args.cores
        t0 = time.time()
        if pdbs:
            print("\nConverting the following to .gro files:")
            for pdb in pdbs:
                print("    "+str(pdb))
            tools.Convert_Structures.pdb2gro(pdbs, cores)
        if gros:
            print("\nConverting the following to .pdb files:")
            for gro in gros:
                print("    "+str(gro))
            tools.Convert_Structures.gro2pdb(gros, cores)
        if trajectories:
            print("\nProcessing the following trajectories:")
            for trj in trajectories:
                print("    "+str(trj))
            if trjconv == 'backbone':
                gromacs_num = 4
            elif trjconv == 'masses':
                gromacs_num = 10
            elif trjconv == 'system':
                gromacs_num = 0
            tools.Convert_Structures.trjconv(
                trajectories, cores=cores, structure_file=structure,
                gromacs_num=gromacs_num)
        t1 = time.time()
        print("\nFinished in "+str(t1-t0)+" seconds\n")
