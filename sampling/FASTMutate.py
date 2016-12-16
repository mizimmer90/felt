# Author: Maxwell Zimmerman <mizimmer@wustl.edu>
# Copyright (c) 2016, Washington University
# All rights reserved.

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from __future__ import absolute_import, print_function, division
import numpy as np
import os
import sys
import time
from .core_sampling import core_sampling

__all__=['FASTMutate']

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

class FASTMutate(core_sampling):
    """FASTMutate

    Samples using populations and PageRanks.

    Parameters
    ----------
    alpha_beta : str, optional, default: 1,1
        The alpha and beta values for ranking. This will control the balance
        between exploration and exploitation. i.e. ri = (alpha*phi)+(beta*psi)
    landscape : str, optional, default: landscape.dat
        The name of the landscape info file. This information will determine
        the objective function.
    statistics : str, optional, default: counts
        The statistical method to use when ranking states. Default is counts,
        although PageRanking methods are available. Choices are:
        1) 'counts'
        2) 'pagerank' : The page rank algorithm applied to state-space.
            PR(i) = (1-d)/N + d*Aij*PR where PR is the page rank, Aij is the
            adjacency matrix, d is the dampening factor, and N is the number
            of nodes. Aij = sum(1/Oj) where Oj is the number of states j
            is observed to transision to. d is taken to be = 0.9.
        3) 'none' : No statistical component

    """
    def __init__(
            self, alpha_beta='1,1', landscape='landscape.dat',
            statistics='counts'):
        self.alpha_beta = self.args.alpha_beta
        self.landscape = self.args.landscape
        self.statistics = self.args.statistics
        core_sampling.__init__(self)

    def print_and_check_inputs(self):
        core_sampling.print_basic_inputs(self)
        print("FAST Alpha,Beta:\t\t"+str(self.alpha_beta))
        print("Landscape filename:\t\t"+str(self.landscape))
        print("FAST-statistics:\t\t"+str(self.statistics))
        print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("~~~~                                                            ~~~~")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
        print
        core_sampling.load_basic_inputs(self)
        return

def load_landscape(filename):
    landscape = np.loadtxt(filename)
    return landscape
