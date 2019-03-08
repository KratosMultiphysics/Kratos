from __future__ import print_function, absolute_import, division
## @module aitken
# This module contains the class Aitken
# Author: Wei He
# Date: Feb. 20, 2017

# Importing the base class
from co_simulation_base_convergence_accelerator import CoSimulationBaseConvergenceAccelerator

# Other imports
import numpy as np
from copy import deepcopy
from collections import deque
from co_simulation_tools import red, classprint

def Create(settings, solvers, level):
    return Aitken(settings, solvers, level)

## Class Aitken.
# This class contains the implementation of Aitken relaxation and helper functions.
# Reference: Ulrich Küttler et al., "Fixed-point fluid–structure interaction solvers with dynamic relaxation"
class Aitken(CoSimulationBaseConvergenceAccelerator):
    ## The constructor.
    # @param init_alpha Initial relaxation factor in the first time step.
    # @param init_alpha_max Maximum relaxation factor for the first iteration in each time step
    def __init__( self, settings, solvers, level ):
        super(Aitken, self).__init__(settings, solvers, level)
        self.R = deque( maxlen = 2 )
        if "init_alpha" in self.settings:
            self.alpha_old = self.settings["init_alpha"]
        else:
            self.alpha_old = 0.1
        if "init_alpha_max" in self.settings:
            self.init_alpha_max = self.settings["init_alpha_max"]
        else:
            self.init_alpha_max = 0.45


    ## ComputeUpdate(r, x)
    # @param r residual r_k
    # @param x solution x_k
    # Computes the approximated update in each iteration.
    def _ComputeUpdate( self, r, x ):
        self.R.appendleft( deepcopy(r) )
        k = len( self.R ) - 1
        ## For the first iteration, do relaxation only
        if k == 0:
            alpha = min( self.alpha_old, self.init_alpha_max )
            if self.echo_level > 3:
                classprint(self.lvl, self._Name(), ": Doing relaxation in the first iteration with initial factor = " + "{0:.1g}".format(alpha))
            return alpha * r
        else:
            r_diff = self.R[0] - self.R[1]
            numerator = np.inner( self.R[1], r_diff )
            denominator = np.inner( r_diff, r_diff )
            alpha = -self.alpha_old * numerator/denominator
            if self.echo_level > 3:
                classprint(self.lvl, self._Name(), ": Doing relaxation with factor = " + "{0:.1g}".format(alpha))
            if alpha > 20:
                alpha = 20
                if self.echo_level > 0:
                    classprint(self.lvl, self._Name(), ": "+ red("WARNING: dynamic relaxation factor reaches upper bound: 20"))
            elif alpha < -2:
                alpha = -2
                if self.echo_level > 0:
                    classprint(self.lvl, self._Name(), ": " + red("WARNING: dynamic relaxation factor reaches lower bound: -2"))
            delta_x = alpha * self.R[0]
        self.alpha_old = alpha

        return delta_x

    def _Name(self):
        return self.__class__.__name__