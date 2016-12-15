# Author: Maxwell I. Zimmerman
# Copyright (c) 2016, Washington University
# Washington University in St. Louis

##############################################################################
# Imports
##############################################################################

import numpy as np
import time

##############################################################################
# Code
##############################################################################

def output_status(comment):
    t0 = time.localtime()
    ftime = "[%d:%d:%d]" % (t0.tm_hour, t0.tm_min, t0.tm_sec)
    print(ftime+" "+str(comment)+"\n")
    return
