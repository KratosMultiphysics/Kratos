# This module contains the class AndersonConvergenceAccelerator
# Author: Andreas Winterstein
# Date: Jul. 2018

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_convergence_accelerator import CoSimulationConvergenceAccelerator

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import classprint, bold, red, magenta, blue

# Other imports
import numpy as np
from copy import deepcopy
from collections import deque

def Create(settings, solver_wrapper):
    cs_tools.SettingsTypeCheck(settings)
    raise NotImplementedError("This class needs some updates see MVQN and Aitken")
    return AndersonConvergenceAccelerator(settings, solver_wrapper)

class AndersonConvergenceAccelerator(CoSimulationConvergenceAccelerator):
    ## The constructor.
    # @param iteration_horizon number of values to be stored of last iterations.
    # @param alpha Relaxation factor for computing the update.
    # @param beta weighting factor of constant relaxation
    # @param p factor for switch between constant relaxation and alternating anderson Gauß-Seidel/Jacobian method
    # p = 1 results in the Anderson acceleration and p -> infinity results in constant relaxation
    def __init__( self, settings, solvers, cosim_solver_details ):
        super().__init__(settings, solvers, cosim_solver_details)
        if "iteration_horizon" in self.settings:
            iteration_horizon = self.settings["iteration_horizon"]
        else:
            iteration_horizon = 20
        if "alpha" in self.settings:
            self.alpha = self.settings["alpha"]
        else:
            self.alpha = 0.1
        if "beta" in self.settings:
            self.beta = self.settings["beta"]
        else:
            self.beta = 0.2
        if "p" in self.settings:
            self.p = self.settings["p"]
        else:
            self.p = 2

        self.V = deque( maxlen = iteration_horizon )
        self.W = deque( maxlen = iteration_horizon )

        self.iteration_counter = 0

    def InitializeSolutionStep(self):
        self.iteration_counter = 0

   ## UpdateSolution(r, x)
    # @param r residual r_k
    # @param x solution x_k
    # Computes the approximated update in each iteration.

    def UpdateSolution(self, r, x):

        self.V.appendleft( deepcopy(r) )
        self.W.appendleft( deepcopy(x) )
        row = len(r)
        col = len( self.V ) - 1
        k = col
        if k == 0:
            ## For the first iteration, do relaxation only
            if self.echo_level > 3:
                classprint(self._ClassName(), "Doing relaxation in the first iteration with factor = ", "{0:.1g}".format(self.alpha))
            return self.alpha * r
        else:
            self.F = np.empty( shape = (col, row) ) # will be transposed later
            self.X = np.empty( shape = (col, row) ) # will be transposed later
            for i in range(0, col):
                self.F[i] = self.V[i] - self.V[i + 1]
                self.X[i] = self.W[i] - self.W[i + 1]
            self.F = self.F.T
            self.X = self.X.T


            #compute Moore-Penrose inverse of F^T F
            A = np.linalg.pinv(self.F.T @ self.F)

            switch = (self.iteration_counter + 1)/self.p

            classprint(magenta(self.iteration_counter))

            if switch.is_integer() == True:
                B = self.beta * np.identity(row) - (self.X + self.beta * self.F) @ A @ self.F.T
                if self.echo_level > 3:
                    classprint(self._ClassName(), blue("Compute B with Anderson"))
            else:
                B = self.alpha * np.identity(row)
                if self.echo_level > 3:
                    classprint(self._ClassName(), red("Constant underrelaxtion"))

            delta_x = B @ r

            self.iteration_counter += 1

            return delta_x

    def FinalizeSolutionStep(self):
        self.V.clear()
        self.W.clear()
