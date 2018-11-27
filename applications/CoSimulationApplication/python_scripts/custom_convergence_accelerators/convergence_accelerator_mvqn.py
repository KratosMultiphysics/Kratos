## @module iqnils
# This module contains the class MVQN
# Author: Wei He
# Date: Feb. 20, 2017

try :
    import numpy as np
except ModuleNotFoundError:
    print(cs_tools.bcolors.FAIL + 'Numpy is not available ! Using python default lists for computation'+ cs_tools.bcolors.ENDC)

from base_co_simulation_classes.co_simulation_base_convergence_accelerator import CoSimulationBaseConvergenceAccelerator
import co_simulation_tools as cs_tools
data_structure = cs_tools.cs_data_structure

from copy import deepcopy
from collections import deque

def Create(settings, solver):
    accelerator = MultiVectorQuasiNewtonAccelerator(settings, solver)
    return accelerator


## Class MVQN.
# This class contains the implementation of the MVQN method and helper functions.
# Reference: A.E.J. Bogaers et al. "Quasi-Newton methods for implicit black-box FSI coupling", Computational methods in applied mechanics and engineering. 279(2014) 113-132.
class MultiVectorQuasiNewtonAccelerator(CoSimulationBaseConvergenceAccelerator):
    ## The constructor.
    # @param horizon Maximum number of vectors to be stored in each time step.
    # @param alpha Relaxation factor for computing the update, when no vectors available.
    def __init__( self, settings, solver ):
        super(MultiVectorQuasiNewtonAccelerator, self).__init__(settings, solver)
        default_settings = self._GetDefaultSettings()
        self.settings.RecursivelyValidateAndAssignDefaults(default_settings)
        self.horizon = self.settings["settings"]["horizon"].GetInt()
        self.R = deque( maxlen = self.horizon )
        self.X = deque( maxlen = self.horizon )
        self.alpha = self.settings["settings"]["alpha"].GetDouble()
        self.J = [] # size will be determined when first time get the input vector
        self.J_hat = []
        self.iteration = 0
        self.data_name = self.settings["data"]["data_name"].GetString()

    ## _GetDefaultSettings :  Function to define the default parameters for this class
    #
    def _GetDefaultSettings(self):
        default_setting = data_structure.Parameters("""
            {
                "type"          : "mvqn",
                "data" : {
                        "solver"   : "",
                        "data_name"     : ""
                },
                "settings":{
                    "horizon" : 15,
                    "alpha" : 0.30
                }
            }
        """)
        return default_setting

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
        print( "Number of new modes: ", col )

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


    ## InitializeSolutionStep : Called once at the beginning of the solution step
    #
    def InitializeSolutionStep(self):
        self.iteration = 0
        if self.J == []:
            return

        row = self.J.shape[0]
        col = self.J.shape[1]
        ## Assign J=J_hat
        for i in range(0, row):
            for j in range(0, col):
                self.J[i][j] = self.J_hat[i][j]

        ## Clear the buffer
        if self.R and self.X:
            self.R.clear()
            self.X.clear()

    ## FinalizeSolutionStep : Called once at the end of the solution step
    #
    def FinalizeSolutionStep(self):
        pass

    ## InitializeNonLinearIteration : Function initializes the non linear iteration (coupling iteration)
    #                               Called at the beginning of the nonlinear iteration (coupling iteration)
    #
    def InitializeNonLinearIteration(self):
        self.input_data_current_iter = cs_tools.GetDataAsList(self.solver, self.data_name)


    ## FinalizeNonLinearIteration : Function finalizes the non linear iteration (coupling iteration)
    #                               Called at the end of the nonlinear iteration (coupling iteration)
    #
    def FinalizeNonLinearIteration(self):
        self.output_data_current_iter = cs_tools.GetDataAsList(self.solver, self.data_name)
        residual = self._CalculateResidual()
        self.update = self._CalculateUpdate(residual, self.output_data_current_iter)
        self._ApplyRelaxationToData()
        self.iteration = self.iteration + 1

    ## PrintInfo : Function to print the information of the convergence accelerator
    #
    def PrintInfo(self):
        print(cs_tools.bcolors.HEADER + "This is an object of Multi-Vector Quasi Newton relaxation accelerator. Horizon alpha is ", self.horizon, ", alpha is : ", self.alpha,""+cs_tools.bcolors.ENDC )

    ## Check : Function to Check the setup of the convergence accelerator
    #
    def Check(self):
        pass

    ## _Name : Function to get the name of the convergence accelerator
    #
    def _Name(self):
        return "mvqn"

    ## _CalculateResidual : Calculates residual of the data specified in the settings
    #                       Numpy can be used in the variants of this class.
    #                       residual = output_data_current_iter - input_data_current_iter
    def _CalculateResidual(self):
        if(self.iteration == 0):
            self.output_data_current_iter
            return

        return self._Difference(self.output_data_current_iter , self.input_data_current_iter)

    ## _ApplyRelaxationToData : updates the data with the update calculated
    #
    def _ApplyRelaxationToData(self):
        updated_data = [ input_data + update  for input_data, update in zip(self.input_data_current_iter, self.update)]
        cs_tools.ApplyUpdateToData(self.solver, self.data_name, updated_data)