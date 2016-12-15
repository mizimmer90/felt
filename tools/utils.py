# Author: Maxwell I. Zimmerman
# Copyright (c) 2016, Washington University
# Washington University in St. Louis

##############################################################################
# Imports
##############################################################################

import itertools
import numpy as np
import subprocess as sp
from multiprocessing import Pool

##############################################################################
# Code
##############################################################################

def convert_to_string(binary):
    return binary.decode('utf-8')

def _run_command(cmd_info):
    cmd,supress = cmd_info
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    output,err = p.communicate()
    p.terminate()
    if convert_to_string(err) != '' and not supress:
        raise command_error(err)
    return output

def run_commands(cmds, supress=False, cores=1):
    cmd_info = list(zip(cmds,itertools.repeat(supress)))
    pool = Pool(processes=cores)
    outputs = pool.map(_run_command, cmd_info)
    pool.terminate()
    return outputs

