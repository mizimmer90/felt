# Author: Maxwell I. Zimmerman
# Copyright (c) 2016, Washington University
# Washington University in St. Louis

##############################################################################
# Imports
##############################################################################

import numpy as np
import os
import time
from ..tools import sim_basics

##############################################################################
# Code
##############################################################################

def _format_time():
    t0 = time.localtime()
    ftime = "[%d:%d:%d]" % (t0.tm_hour, t0.tm_min, t0.tm_sec)
    return ftime

def output_status(comment):
    ftime = _format_time()
    print(ftime+" "+str(comment)+"\n")
    return

def output_mutation(res_num,prev_res,new_res):
    ftime = _format_time()
    print(
        ftime+" mutating residue number %d from %s to %s" \
        % (res_num, prev_res, new_res))
    return

def save_fixers_as_pdbs(fixers, output_names=None):
    if output_names is None:
        output_names = np.array(
            ['State%d.pdb' % num for num in range(len(fixers))], dtype=str)
    for num in range(len(fixers)):
        pdb = sim_basics.pdb_from_fixer(fixers[num])
        pdb.save_pdb(output_names[num])
    return
