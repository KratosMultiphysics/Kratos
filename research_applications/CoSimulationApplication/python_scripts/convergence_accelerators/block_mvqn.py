## @module block_mvqn
# This module contains the class MVQNConvergenceAccelerator
# Author: Tiba Azzeddine
# Date: Nov. 06, 2023

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_convergence_accelerator import CoSimulationConvergenceAccelerator

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# Other imports
import numpy as np
from copy import deepcopy
from collections import deque


def Create(settings):
    cs_tools.SettingsTypeCheck(settings)
    return BLOCKMVQNConvergenceAccelerator(settings)

## Class BLOCKMVQNConvergenceAccelerator.
# This class contains the implementation of the Block MVQN method and helper functions.
# Reference: A.E.J. Bogaers et al. "Quasi-Newton methods for implicit black-box FSI coupling", Computational methods in applied mechanics and engineering. 279(2014) 113-132.
class BLOCKMVQNConvergenceAccelerator(CoSimulationConvergenceAccelerator):
    ## The constructor.
    # @param horizon Maximum number of vectors to be stored in each time step.
    # @param alpha Relaxation factor for computing the update, when no vectors available.
    def __init__( self, settings):
        super().__init__(settings)

        horizon = self.settings["horizon"].GetInt()
        self.alpha = self.settings["alpha"].GetDouble()
        self.epsilon = self.settings["epsilon"].GetDouble()

        self.X_tilde = {}
        self.X = {}
        self.J = {}
        self.J_hat = {}
        self.coupl_data_names = {}

        for solver_data in settings["solver_sequence"].values():
            self.X_tilde[solver_data["data_name"].GetString()] = deque( maxlen = horizon )
            self.X[solver_data["data_name"].GetString()] = deque( maxlen = horizon )
            self.J[solver_data["data_name"].GetString()] = None # size will be determined when first time get the input vector
            self.J_hat[solver_data["data_name"].GetString()] = None
            self.coupl_data_names[solver_data["data_name"].GetString()] = solver_data["coupled_data_name"].GetString()

    ## UpdateSolution(r, x, y, data_name, yResidual)
    # @param r residual r_k
    # @param x solution x_k
    # @param y (coupled solver) solution y_k
    # @param data_name coupling variable
    # @param yResidual (coupled solver) residual yResidual
    # Computes the approximated update in each iteration.
    def UpdateSolution( self, r, x, y, data_name, yResidual,):

        coupled_data_name = self.coupl_data_names[data_name]
        self.X_tilde[data_name].appendleft( deepcopy(r + x) )
        self.X[coupled_data_name].appendleft( deepcopy(y) )

        col = len(self.X[coupled_data_name]) - 1
        row = len(r)
        rowY = len(y)
        k = col
        if self.echo_level > 3:
            cs_tools.cs_print_info(self._ClassName(), "Number of new modes: ", col )

        ## For the first iteration
        if k == 0:
            if self.J[data_name] is None or self.J[coupled_data_name] is None:
                ## Zero initial Jacobians
                self.J_hat[data_name] = np.zeros( (row, rowY) )
                self.J_hat[coupled_data_name] = np.zeros( (rowY, row) )
                self.J[coupled_data_name] = np.zeros( (rowY, row) )
                self.J[data_name] = np.zeros( (row, rowY) )
                return self.alpha * r # Initial acceleration to be a constant relaxation
            else:

                blockJacobian = (np.eye(row) - self.J[data_name] @ self.J[coupled_data_name])
                b = r - self.J[data_name] @ yResidual
                return np.linalg.solve( blockJacobian, b )

        ## Construct matrix W(differences of intermediate solutions y)
        W = np.empty( shape = (col, rowY) ) # will be transposed later
        for i in range(0, col):
            W[i] = self.X[coupled_data_name][i] - self.X[coupled_data_name][i + 1]
        W = W.T

        ## Construct matrix W(differences of intermediate solutions x~)
        V = np.empty( shape = (col, row) ) # will be transposed later
        for i in range(0, col):
            V[i] = self.X_tilde[data_name][i] - self.X_tilde[data_name][i + 1]
        V = V.T

        self.J_hat[data_name] = self.J[data_name] + (V - self.J[data_name] @ W) @ (np.linalg.pinv(W, rcond=self.epsilon))

        blockJacobian = (np.eye(row) - self.J_hat[data_name] @ self.J_hat[coupled_data_name])
        b = r - self.J_hat[data_name] @ yResidual

        return np.linalg.solve( blockJacobian, b )

    def FinalizeSolutionStep( self ):

        ## Assign J=J_hat
        for data_name in self.J:
            self.J[data_name] = self.J_hat[data_name].copy()

            ## Clear the buffer
            self.X_tilde[data_name].clear()
            self.X[data_name].clear()

        if self.echo_level > 3:
            cs_tools.cs_print_info(self._ClassName(), "Jacobian matrix updated!")

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "horizon" : 15,
            "alpha"   : 1.0,
            "epsilon" : 1e-9,
            "solver_sequence" : []
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
