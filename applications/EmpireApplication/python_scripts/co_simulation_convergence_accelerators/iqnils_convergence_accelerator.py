from __future__ import print_function, absolute_import, division
## @module iqnils
# This module contains the class IQNILS
# Author: Wei He
# Date: Feb. 20, 2017

# Importing the base class
from co_simulation_base_convergence_accelerator import CoSimulationBaseConvergenceAccelerator

# Other imports
import numpy as np
from copy import deepcopy
from collections import deque
from co_simulation_tools import classprint

def Create(settings, solvers, level):
    return IQNILS(settings, solvers, level)

## Class IQNILS.
# This class contains the implementation of the IQN-ILS method and helper functions.
# Reference: Joris Degroote, PhD thesis "Development of algorithms for the partitioned simulation of strongly coupled fluid-structure interaction problems", 84-91.
class IQNILS(CoSimulationBaseConvergenceAccelerator):
    ## The constructor.
    # @param iteration_horizon Maximum number of vectors to be stored in each time step.
    # @param timestep_horizon Maximum number of time steps of which the vectors are used.
    # @param alpha Relaxation factor for computing the update, when no vectors available.
    def __init__( self, settings, solvers, level ):
        super(IQNILS, self).__init__(settings, solvers, level)
        if "iteration_horizon" in self.settings:
            iteration_horizon = self.settings["iteration_horizon"]
        else:
            iteration_horizon = 20
        if "timestep_horizon" in self.settings:
            timestep_horizon = self.settings["timestep_horizon"]
        else:
            timestep_horizon = 1
        if "alpha" in self.settings:
            self.alpha = self.settings["alpha"]
        else:
            self.alpha = 0.125

        self.R = deque( maxlen = iteration_horizon )
        self.X = deque( maxlen = iteration_horizon )
        self.q = timestep_horizon - 1
        self.v_old_matrices = deque( maxlen = self.q )
        self.w_old_matrices = deque( maxlen = self.q )
        self.V_new = []
        self.W_new = []
        self.V_old = []
        self.W_old = []

    ## ComputeUpdate(r, x)
    # @param r residual r_k
    # @param x solution x_k
    # Computes the approximated update in each iteration.
    def _ComputeUpdate( self, r, x ):
        self.R.appendleft( deepcopy(r) )
        self.X.appendleft(    x + r    )  # r = x~ - x
        row = len(r)
        col = len( self.R ) - 1
        k = col
        num_old_matrices = len( self.v_old_matrices )

        if self.V_old == [] and self.W_old == []: # No previous vectors to reuse
            if k == 0:
                ## For the first iteration in the first time step, do relaxation only
                if self.echo_level > 3:
                    classprint(self.lvl, self._Name(), "Doing relaxation in the first iteration with factor = ", "{0:.1g}".format(self.alpha))
                return self.alpha * r
            else:
                if self.echo_level > 3:
                    classprint(self.lvl, self._Name(), "Doing multi-vector extrapolation")
                    classprint(self.lvl, self._Name(), "Number of new modes: ", col)
                self.V_new = np.empty( shape = (col, row) ) # will be transposed later
                for i in range(0, col):
                    self.V_new[i] = self.R[i] - self.R[i + 1]
                self.V_new = self.V_new.T
                V = self.V_new

                ## Check the dimension of the newly constructed matrix
                if ( V.shape[0] < V.shape[1] ) and self.echo_level > 0:
                    classprint(self.lvl, self._Name(), ": "+ red("WARNING: column number larger than row number!"))

                ## Construct matrix W(differences of predictions)
                self.W_new = np.empty( shape = (col, row) ) # will be transposed later
                for i in range(0, col):
                    self.W_new[i] = self.X[i] - self.X[i + 1]
                self.W_new = self.W_new.T
                W = self.W_new

                ## Solve least-squares problem
                delta_r = -self.R[0]
                c = np.linalg.lstsq(V, delta_r)[0]

                ## Compute the update
                delta_x = np.dot(W, c) - delta_r

                return delta_x
        else:  # previous vectors can be reused
            if k == 0: # first iteration
                if self.echo_level > 3:
                    classprint(self.lvl, self._Name(), "Using matrices from previous time steps")
                    classprint(self.lvl, self._Name(), "Number of previous matrices: ", num_old_matrices)
                V = self.V_old
                W = self.W_old
                ## Solve least-squares problem
                delta_r = -self.R[0]
                c = np.linalg.lstsq(V, delta_r)[0]

                ## Compute the update
                delta_x = np.dot(W, c) - delta_r
                return delta_x
            else:
                ## For other iterations, construct new V and W matrices and combine them with old ones
                if self.echo_level > 3:
                    classprint(self.lvl, self._Name(), "Doing multi-vector extrapolation")
                    classprint(self.lvl, self._Name(), "Number of new modes: ", col)
                    classprint(self.lvl, self._Name(), "Number of previous matrices: ", num_old_matrices)
                ## Construct matrix V (differences of residuals)
                self.V_new = np.empty( shape = (col, row) ) # will be transposed later
                for i in range(0, col):
                    self.V_new[i] = self.R[i] - self.R[i + 1]
                self.V_new = self.V_new.T
                V = np.hstack( (self.V_new, self.V_old) )
                ## Check the dimension of the newly constructed matrix
                if ( V.shape[0] < V.shape[1] ) and self.echo_level > 0:
                    classprint(self.lvl, self._Name(), ": "+ red("WARNING: column number larger than row number!"))

                ## Construct matrix W(differences of predictions)
                self.W_new = np.empty( shape = (col, row) ) # will be transposed later
                for i in range(0, col):
                    self.W_new[i] = self.X[i] - self.X[i + 1]
                self.W_new = self.W_new.T
                W = np.hstack( (self.W_new, self.W_old) )

                ## Solve least-squares problem
                delta_r = -self.R[0]
                c = np.linalg.lstsq(V, delta_r)[0]

                ## Compute the update
                delta_x = np.dot(W, c) - delta_r

                return delta_x

    ## FinalizeSolutionStep()
    # Finalizes the current time step and initializes the next time step.
    def FinalizeSolutionStep( self ):
        if self.V_new != [] and self.W_new != []:
            self.v_old_matrices.appendleft( self.V_new )
            self.w_old_matrices.appendleft( self.W_new )
        if self.v_old_matrices and self.w_old_matrices:
            self.V_old = np.concatenate( self.v_old_matrices, 1 )
            self.W_old = np.concatenate( self.w_old_matrices, 1 )
        ## Clear the buffer
        if self.R and self.X:
            if self.echo_level > 3:
                classprint(self.lvl, self._Name(), "Cleaning")
            self.R.clear()
            self.X.clear()
        self.V_new = []
        self.W_new = []

    def _Name(self):
        return self.__class__.__name__