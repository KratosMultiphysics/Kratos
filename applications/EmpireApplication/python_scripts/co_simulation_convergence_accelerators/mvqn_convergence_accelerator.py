from __future__ import print_function, absolute_import, division
## @module iqnils
# This module contains the class MVQN
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
    return MVQN(settings, solvers, level)

## Class MVQN.
# This class contains the implementation of the MVQN method and helper functions.
# Reference: A.E.J. Bogaers et al. "Quasi-Newton methods for implicit black-box FSI coupling", Computational methods in applied mechanics and engineering. 279(2014) 113-132.
class MVQN(CoSimulationBaseConvergenceAccelerator):
    ## The constructor.
    # @param horizon Maximum number of vectors to be stored in each time step.
    # @param alpha Relaxation factor for computing the update, when no vectors available.
    def __init__( self, settings, solvers, level ):
        super(MVQN, self).__init__(settings, solvers, level)
        if "horizon" in self.settings:
            horizon = self.settings["horizon"]
        else:
            horizon = 15
        if "alpha" in self.settings:
            self.alpha = self.settings["alpha"]
        else:
            self.alpha = 0.125
        self.R = deque( maxlen = horizon )
        self.X = deque( maxlen = horizon )
        self.J = [] # size will be determined when first time get the input vector
        self.J_hat = []

    ## ComputeUpdate(r, x)
    # @param r residual r_k
    # @param x solution x_k
    # Computes the approximated update in each iteration.
    def _ComputeUpdate( self, r, x ):
        self.R.appendleft( deepcopy(r) )
        self.X.appendleft( deepcopy(x) )
        col = len(self.R) - 1
        row = len(r)
        k = col
        if self.echo_level > 3:
            classprint(self.lvl, self._Name(), "Number of new modes: ", col )

        ## For the first iteration
        if k == 0:
            if self.J == []:
                return self.alpha * r  # if no Jacobian, do relaxation
            else:
                return np.linalg.solve( self.J, -r ) # use the Jacobian from previous step

        ## Let the initial Jacobian correspond to a constant relaxation
        if self.J == []:
            self.J = - np.identity( row ) / self.alpha # correspongding to constant relaxation

        ## Construct matrix V (differences of residuals)
        V = np.empty( shape = (col, row) ) # will be transposed later
        for i in range(0, col):
            V[i] = self.R[i] - self.R[i + 1]
        V = V.T

        ## Construct matrix W(differences of intermediate solutions x)
        W = np.empty( shape = (col, row) ) # will be transposed later
        for i in range(0, col):
            W[i] = self.X[i] - self.X[i + 1]
        W = W.T

        ## Solve least norm problem
        rhs = V - np.dot(self.J, W)
        b = np.identity( row )
        W_right_inverse = np.linalg.lstsq(W, b)[0]
        J_tilde = np.dot(rhs, W_right_inverse)
        self.J_hat = self.J + J_tilde
        delta_r = -self.R[0]
        delta_x = np.linalg.solve(self.J_hat, delta_r)

        return delta_x

    ## FinalizeSolutionStep()
    # Finalizes the current time step and initializes the next time step.
    def FinalizeSolutionStep( self ):
        if self.J == []:
            return

        row = self.J.shape[0]
        col = self.J.shape[1]
        ## Assign J=J_hat
        for i in range(0, row):
            for j in range(0, col):
                self.J[i][j] = self.J_hat[i][j]
        if self.echo_level > 3:
            classprint(self.lvl, self._Name(), "Jacobian matrix updated!")
        ## Clear the buffer
        if self.R and self.X:
            self.R.clear()
            self.X.clear()

    def _Name(self):
        return self.__class__.__name__