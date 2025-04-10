## @module iqnils
# This module contains the class IQNMVJExplicitConvergenceAccelerator
# Author: Wei He
# Date: Feb. 20, 2017
# Revised: Susanna Baars
# Date: April 10, 2025

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
import typing

def Create(settings):
    cs_tools.SettingsTypeCheck(settings)
    return IQNIMVJExplicitConvergenceAccelerator(settings)

## Class IQNMVJExplicitConvergenceAccelerator.
# This class contains the implementation of the explicit IQN-MVJ method and helper functions.
# Reference: F. Lindner et al. "A comparison of various quasi-Newton schemes for partitioned fluid-structure interaction", 6th International Conference on Computational Methods for Coupled Problems in Science and Engineering; 2015. p. 477–88.
class IQNIMVJExplicitConvergenceAccelerator(CoSimulationConvergenceAccelerator):
    ## The constructor.
    # @param iteration_horizon Maximum number of vectors to be stored in each time step.
    # @param alpha Relaxation factor for computing the update, when no vectors available.
    def __init__( self, settings):
        super().__init__(settings)

        iteration_horizon = self.settings["iteration_horizon"].GetInt()
        self.alpha = self.settings["alpha"].GetDouble()

        self.R = deque( maxlen = iteration_horizon )
        self.X = deque( maxlen = iteration_horizon )
        self.J_inv: typing.Optional[np.ndarray] = None # size will be determined when first time get the input vector
        self.J_inv_hat: typing.Optional[np.ndarray] = None

    @staticmethod
    def _construct_difference_matrices(R, X, col):
        V = np.array([R[i + 1] - R[0] for i in range(col)]).T
        W = np.array([X[i + 1] - X[0] for i in range(col)]).T
        
        return V, W

    ## UpdateSolution(r, x)
    # @param r residual r_k
    # @param x input x_k
    # Computes the approximated update in each iteration.
    def UpdateSolution( self, r, x ):
        self.R.appendleft(deepcopy(r))
        self.X.appendleft(x + r)  # r = x~ - x
        col = len(self.R) - 1
        row = len(r)
        k = col
        delta_r = - r

        ## For the first iteration
        if k == 0:
            if self.J_inv is None:
                if self.echo_level > 3:
                    cs_tools.cs_print_info(self._ClassName(), "Doing relaxation in the first iteration with factor = ", "{0:.1g}".format(self.alpha))

                return self.alpha * r  # if no inverse Jacobian, do relaxation
            else:
                if self.echo_level > 3:
                    cs_tools.cs_print_info(self._ClassName(), "Applying the explicit IQN-IMVJ method.")

                delta_x = self.J_inv @ delta_r - delta_r  # use the inverse Jacobian from previous step

                return delta_x

        if self.J_inv is None:
            self.J_inv = np.zeros([row, row]) # Initialize the inverse Jacobian as zero

        # Prepare inputs and outputs
        V, W  = self._construct_difference_matrices(self.R, self.X, col)
        W_diff = W - self.J_inv @ V

        # Update the Jacobian
        V_right_inverse = np.linalg.lstsq(V, np.identity(row), rcond=-1)[0]
        J_update = W_diff @ V_right_inverse

        K = np.dot(V.T, V)
        try:
            Z = np.linalg.solve(K, V.T)
        except:
            Z = np.linalg.lstsq(K, V.T)[0]

        J_update = W_diff @ Z

        self.J_inv_hat = self.J_inv + J_update

        # Compute update
        delta_x = self.J_inv_hat @ delta_r - delta_r

        return delta_x

    ## FinalizeSolutionStep()
    # Finalizes the current time step and initializes the next time step.
    def FinalizeSolutionStep( self ):
        if self.J_inv_hat is not None:
            ## Assign J_inv = J_inv_hat
            self.J_inv = np.copy(self.J_inv_hat)
        if self.echo_level > 3:
            cs_tools.cs_print_info(self._ClassName(), "Jacobian matrix updated!")
        ## Clear the buffer
        self.R.clear()
        self.X.clear()

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "iteration_horizon" : 15,
            "alpha"   : 0.125
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
