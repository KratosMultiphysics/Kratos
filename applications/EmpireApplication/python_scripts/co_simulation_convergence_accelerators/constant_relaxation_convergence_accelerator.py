from __future__ import print_function, absolute_import, division

# Importing the base class
from co_simulation_base_convergence_accelerator import CoSimulationBaseConvergenceAccelerator

# Other imports
import numpy as np
from copy import deepcopy
from collections import deque
from co_simulation_tools import classprint

def Create(settings, solvers, cosim_solver_details, level):
    return ConstantRelaxation(settings, solvers, cosim_solver_details, level)

class ConstantRelaxation(CoSimulationBaseConvergenceAccelerator):
    ## The constructor.
    # @param alpha relaxation factor.
    def __init__( self, settings, solvers, cosim_solver_details, level ):
        super(ConstantRelaxation, self).__init__(settings, solvers, cosim_solver_details, level)
        if "alpha" in self.settings:
            self.alpha = self.settings["alpha"]
        else:
            self.alpha = 0.125

    ## ComputeUpdate(r, x)
    # @param r residual r_k
    # @param x solution x_k
    # Computes the approximated update in each iteration.
    def _ComputeUpdate( self, r, x ):
        if self.echo_level > 3:
            classprint(self.lvl, self._Name(), "Doing relaxation with factor = ", "{0:.1g}".format(self.alpha))
        delta_x = self.alpha * r
        return delta_x

    def _Name(self):
        return self.__class__.__name__