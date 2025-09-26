## @module iqnils
# This module contains the class IQNMVJImplicitConvergenceAccelerator
# Author: Susanna Baars
# Date: April 10, 2025

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_convergence_accelerator import CoSimulationConvergenceAccelerator

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import cs_print_info, cs_print_warning, SettingsTypeCheck
import KratosMultiphysics.CoSimulationApplication.colors as colors

# Other imports
import numpy as np
from copy import deepcopy
from collections import deque

def Create(settings):
    SettingsTypeCheck(settings)
    return IQNIMVJImplicitConvergenceAccelerator(settings)

class _RidgeRegression:
    def __init__(self, regularization_param=1e-14):
        self.regularization_param = regularization_param

    def fit(self, X, y):
        self.X = np.copy(X)
        self.y = np.copy(y)

        K = np.dot(X, X.T)
        n_samples = K.shape[0]
        K.flat[:: n_samples + 1] += self.regularization_param
        try:
            self.Z = np.linalg.solve(K, X).T
        except:
            self.Z = np.linalg.lstsq(K, X)[0].T

    def predict(self, X_test):
        alpha = np.dot(X_test, self.Z)
        y_test = np.dot(alpha, self.y)

        return y_test

## Class IQNMVJImplicitConvergenceAccelerator.
# This class contains the implementation of the IQN-ILS method and helper functions.
# Reference: T. Spenke et al. "A multi-vector interface quasi-Newton method with linear complexity for partitioned fluid–structure interaction", Computer Methods in Applied Mechanics and Engineering, Volume 361; 2020
class IQNIMVJImplicitConvergenceAccelerator(CoSimulationConvergenceAccelerator):
    # @param iteration_horizon Maximum number of vectors to be stored in each time step.
    # @param timestep_horizon Maximum number of time steps of which the vectors are used.
    # @param alpha Relaxation factor for computing the update, when no vectors available.
    def __init__( self, settings):
        super().__init__(settings)

        iteration_horizon = self.settings["iteration_horizon"].GetInt()
        timestep_horizon = self.settings["timestep_horizon"].GetInt()
        self.alpha = self.settings["alpha"].GetDouble()
        self.regularization_param = self.settings["regularization_param"].GetDouble()

        self.R = deque(maxlen = iteration_horizon)
        self.X = deque(maxlen = iteration_horizon)
        self.q = timestep_horizon
        
        self.v_matrices = deque(maxlen = self.q)
        self.w_matrices = deque(maxlen = self.q)
        
        self.rr = []
        
    @staticmethod
    def _make_multifidelity_prediction(rr, x):
        if not rr:
            return np.zeros_like(x)

        return np.sum([_rr.predict(x).reshape(np.shape(x)) for _rr in rr], axis=0)
        
    def _add_fidelity_level(self, V, W):
        W_diff = W - self._make_multifidelity_prediction(self.rr, V.T).T
        rr = _RidgeRegression(regularization_param=self.regularization_param)
        rr.fit(V.T, W_diff.T)
        self.rr.append(rr)
        
    def _compute_update(self, r, k):
        delta_r = -r
                    
        if k == 0 and len(self.v_matrices) == self.v_matrices.maxlen and len(self.w_matrices) == self.w_matrices.maxlen:
            self.rr = []
            for i in range(len(self.v_matrices)):
                V = self.v_matrices[len(self.v_matrices) - 1 - i]
                W = self.w_matrices[len(self.v_matrices) - 1 - i]
                if not len(W) == 0 and not len(V) == 0:
                    self._add_fidelity_level(V, W)

        if k > 1:
            self.rr = self.rr[:-1]

        if k > 0:
            W = self.w_matrices[0]
            V = self.v_matrices[0]
            self._add_fidelity_level(V, W)
        
        delta_x = self._make_multifidelity_prediction(self.rr, delta_r.reshape(1, -1))[0] - delta_r

        return delta_x
    
    @staticmethod
    def _construct_difference_matrices(R, X, col):
        V = np.array([R[i + 1] - R[0] for i in range(col)]).T
        W = np.array([X[i + 1] - X[0] for i in range(col)]).T
        
        return V, W
    
    def _update_stored_matrices(self, V, W, k):
        if k > 0:
            self.v_matrices.popleft()
            self.w_matrices.popleft()

        self.v_matrices.appendleft(V)
        self.w_matrices.appendleft(W)

    ## UpdateSolution(r, x)
    # @param r residual r_k
    # @param x input x_k
    # Computes the approximated update in each iteration.
    def UpdateSolution(self, r, x):
        self.R.appendleft(deepcopy(r))
        self.X.appendleft(x + r)  # r = x_tilde - x
        col = len(self.R) - 1
        k = col

        V, W  = self._construct_difference_matrices(self.R, self.X, col)
        self._update_stored_matrices(V, W, k)

        if k == 0 and all(arr.size == 0 for arr in self.v_matrices) and all(arr.size == 0 for arr in self.w_matrices):
            if self.echo_level > 3:
                cs_print_info(self._ClassName(), "Doing relaxation in the first iteration with factor = ", "{0:.1g}".format(self.alpha))

            return self.alpha * r
        else:
            if self.echo_level > 3:
                cs_print_info(self._ClassName(), "Applying the implicit IQN-IMVJ method.")
            
            delta_x = self._compute_update(r, k)
            
            return delta_x

    # Finalizes the current time step and initializes the next time step.
    def FinalizeSolutionStep( self ):
        ## Clear the buffer
        if self.R and self.X:
            if self.echo_level > 3:
                cs_print_info(self._ClassName(), "Cleaning")
            self.R.clear()
            self.X.clear()

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "iteration_horizon"               : 20,
            "timestep_horizon"                : 50,
            "alpha"                           : 0.125,
            "regularization_param"            : 0.0
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
