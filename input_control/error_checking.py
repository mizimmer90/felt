# Author: Maxwell Zimmerman <mizimmer@wustl.edu>
# Copyright (c) 2016, Washington University
# All rights reserved.

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import numpy as np
import os
import sys

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

allowable_exts = ['gro','pdb']

def check_filenames(filenames):
    if type(filenames) == str:
        filenames = np.loadtxt(filenames)
    for filename in filenames:
        print(filename)    
