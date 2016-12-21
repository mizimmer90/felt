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

def output_status(comment):
    t0 = time.localtime()
    ftime = "[%d:%d:%d]" % (t0.tm_hour, t0.tm_min, t0.tm_sec)
    print(ftime+" "+str(comment)+"\n")
    return

def save_fixers_as_pdbs(fixers, output_names=None):
    if output_names is None:
        output_names = np.array(
            ['State%d.pdb' % num for num in range(len(fixers))], dtype=str)
    for num in range(len(fixers)):
        pdb = sim_basics.pdb_from_fixer(fixers[num])
        pdb.save_pdb(output_names[num])
    return
