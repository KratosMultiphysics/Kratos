## @module ibqnls
# This module contains the class IBQNLSConvergenceAccelerator
# Author: Tiba Azzeddine
# Date: Nov. 13, 2023

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
import scipy as sp


def Create(settings):
    cs_tools.SettingsTypeCheck(settings)
    return IBQNLSConvergenceAccelerator(settings)

## Class IBQNLSConvergenceAccelerator.
# This class contains the implementation of the IBQNLS method and helper functions.
# Reference: Vierendeels J, Lanoye L, Degroote J, Verdonck PR - Implicit coupling of partitioned fluid–structure interaction
#  problems with reduced order models. Comput Struct (2007)
class IBQNLSConvergenceAccelerator(CoSimulationConvergenceAccelerator):
    def __init__( self, settings):
        super().__init__(settings)

        horizon = self.settings["iteration_horizon"].GetInt()
        timestep_horizon = self.settings["timestep_horizon"].GetInt()
        self.alpha = self.settings["alpha"].GetDouble()
        self.gmres_rel_tol = self.settings["gmres_rel_tol"].GetDouble()
        self.gmres_abs_tol = self.settings["gmres_abs_tol"].GetDouble()
        self.q = timestep_horizon - 1

        self.X_tilde = {}
        self.prevX_tilde = None
        self.X = {}
        self.v_old_matrices = {}
        self.w_old_matrices = {}
        self.coupl_data_names = {}
        self.V_old = {}
        self.W_old = {}
        self.W_new = {}
        self.V_new = {}
        self.previousV = None
        self.previousQ = None
        self.previousR = None

        for solverData in settings["solver_sequence"].values():
            self.X_tilde[solverData["data_name"].GetString()] = deque( maxlen = horizon )
            self.X[solverData["data_name"].GetString()] = deque( maxlen = horizon )
            self.v_old_matrices[solverData["data_name"].GetString()] = deque( maxlen = self.q )
            self.w_old_matrices[solverData["data_name"].GetString()] = deque( maxlen = self.q )
            self.V_old[solverData["data_name"].GetString()] = []
            self.W_old[solverData["data_name"].GetString()] = []
            self.V_new[solverData["data_name"].GetString()] = []
            self.W_new[solverData["data_name"].GetString()] = []
            self.coupl_data_names[solverData["data_name"].GetString()] = solverData["coupled_data_name"].GetString()

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

        isFirstDt = (self.V_old[data_name] == [] and self.W_old [data_name] == [])
        isFirstiter = (k == 0)
        # ================== First time step ==================================================
        if isFirstDt:
            # ------------------ First iteration (First time step) ---------------------
            if isFirstiter:
                return self.alpha * r # Initial acceleration to be a constant relaxation

        # ------------------ First iteration (Other time steps) ------------------------
        if isFirstiter:
            V = self.V_old[data_name]
            W = self.W_old[data_name]

        else:
            ## Construct matrix W(differences of intermediate solutions x)
            self.W_new[data_name] = np.empty( shape = (col, rowY) ) # will be transposed later
            for i in range(0, col):
                self.W_new[data_name][i] = self.X[coupled_data_name][i] - self.X[coupled_data_name][i + 1]
            self.W_new[data_name] = self.W_new[data_name].T
            W = self._augmentedMatrix(self.W_new[data_name], self.W_old[data_name], isFirstDt)

            ## Construct matrix W(differences of intermediate solutions y~)
            self.V_new[data_name] = np.empty( shape = (col, row) ) # will be transposed later
            for i in range(0, col):
                self.V_new[data_name][i] = self.X_tilde[data_name][i] - self.X_tilde[data_name][i + 1]
            self.V_new[data_name] = self.V_new[data_name].T
            V = self._augmentedMatrix(self.V_new[data_name], self.V_old[data_name], isFirstDt)

        Q, R = np.linalg.qr(W)
        ##TODO QR Filtering
        b = r - self.pinvProduct(V, Q, R, yResidual)

        if self.previousQ is not None:
            ## Retrieving previous data for the coupled data jacobian approximation ----------------------------------------
            previousQ = self.previousQ
            previousR = self.previousR
            if self.previousV is not None:
                previousV = self.previousV
            else:
                if isFirstDt:
                    previousV = self.V_new[coupled_data_name].copy()
                else:
                    previousV = np.hstack((self.V_new[coupled_data_name].copy(), self.V_old[coupled_data_name].copy()))

            ## Matrix-free implementation of the linear solver (and the pseudoinverse) -------------------------------------
            block_oper = lambda vec: vec - self.pinvProduct(V, Q, R, self.pinvProduct(previousV, previousQ, previousR, vec))
            block_x = sp.sparse.linalg.LinearOperator((row, row), block_oper)
            delta_x, _ = sp.sparse.linalg.gmres( block_x, b, atol=self.gmres_abs_tol, tol=self.gmres_rel_tol )
        else:
            ## Using J = 0 if a previous approximate Jacobian is not available
            delta_x = b

        ## Saving data for the approximate jacobian for the next iteration
        self.previousQ = Q.copy()
        self.previousR = R.copy()
        if isFirstiter: # V only needs to be saved completely on the first iteration
            self.previousV = V.copy()
        else: # Otherwise it can be retrieved from V_new and V_old and memory is freed
            self.previousV = None

        return delta_x

    def _augmentedMatrix(self, mat, oldMat, isFirstDt):
        if isFirstDt:
            return mat.copy()
        else:
            return np.hstack( (mat, oldMat) )

    def pinvProduct(self, LHS, Q, R, x):
        rhs = Q.T @ x
        return LHS @ sp.linalg.solve_triangular(R, rhs)

    def FinalizeSolutionStep( self ):

        for data_name in self.W_new:

            if data_name == list(self.W_new.keys())[-1]: ## Check if last solver in the sequence
                # Saving the V matrix for the next (first) iteration to recover the approximate jacobian
                if self.V_old[data_name] != []:
                    self.previousV = np.hstack((self.V_new[data_name].copy(), self.V_old[data_name].copy()))
                else:
                    self.previousV = self.V_new[data_name].copy()

            if self.V_new[data_name] != [] and self.W_new[data_name] != []:
                self.v_old_matrices[data_name].appendleft( self.V_new[data_name] )
                self.w_old_matrices[data_name].appendleft( self.W_new[data_name] )

            if self.w_old_matrices[data_name] and self.v_old_matrices[data_name]:
                self.V_old[data_name] = np.concatenate( self.v_old_matrices[data_name], 1 )
                self.W_old[data_name] = np.concatenate( self.w_old_matrices[data_name], 1 )

            ## Clear the buffer
            self.X_tilde[data_name].clear()
            self.X[data_name].clear()

        for data_name in self.W_new:
            self.W_new[data_name] = []
            self.V_new[data_name] = []

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "iteration_horizon" : 15,
            "timestep_horizon"  : 3,
            "alpha"   : 1.0,
            "gmres_rel_tol" : 1e-5,
            "gmres_abs_tol" : 1e-14,
            "solver_sequence" : []
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
