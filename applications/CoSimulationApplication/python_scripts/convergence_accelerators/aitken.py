## @module aitken
# This module contains the class AitkenConvergenceAccelerator
# Author: Wei He
# Date: Feb. 20, 2017

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
    return AitkenConvergenceAccelerator(settings)

## Class AitkenConvergenceAccelerator.
# This class contains the implementation of Aitken relaxation and helper functions.
# Reference: Ulrich Küttler et al., "Fixed-point fluid–structure interaction solvers with dynamic relaxation"
class AitkenConvergenceAccelerator(CoSimulationConvergenceAccelerator):
    ## The constructor.
    # @param init_alpha Initial relaxation factor in the first time step.
    # @param init_alpha_max Maximum relaxation factor for the first iteration in each time step
    def __init__( self, settings):
        super().__init__(settings)
        self.R = deque( maxlen = 2 )

        self.alpha_old      = self.settings["init_alpha"].GetDouble()
        self.init_alpha_max = self.settings["init_alpha_max"].GetDouble()
        self.alpha_max      = self.settings["alpha_max"].GetDouble()
        self.alpha_min      = self.settings["alpha_min"].GetDouble()

    def InitializeSolutionStep(self):
        self.initial_iteration = True

    ## UpdateSolution(r, x)
    # @param r residual r_k
    # @param x solution x_k
    # Computes the approximated update in each iteration.
    def UpdateSolution( self, r, x ):
        self.R.appendleft( deepcopy(r) )

        ## For the first iteration, do relaxation only
        if self.initial_iteration:
            self.initial_iteration = False
            alpha = min( self.alpha_old, self.init_alpha_max )
            if self.echo_level > 3:
                cs_tools.cs_print_info(self._ClassName(), ": Doing relaxation in the first iteration with initial factor = {}".format(alpha))
            return alpha * r

        else:
            r_diff = self.R[0] - self.R[1]
            numerator = np.inner( self.R[1], r_diff )
            denominator = np.inner( r_diff, r_diff )
            alpha = -self.alpha_old * numerator/denominator
            if self.echo_level > 3:
                cs_tools.cs_print_info(self._ClassName(), ": Doing relaxation with factor = {}".format(alpha))
            if alpha > self.alpha_max:
                alpha = self.alpha_max
                if self.echo_level > 0:
                    cs_tools.cs_print_warning(self._ClassName(), "dynamic relaxation factor reaches upper bound: {}".format(self.alpha_max))
            elif alpha < self.alpha_min:
                alpha = self.alpha_min
                if self.echo_level > 0:
                    cs_tools.cs_print_warning(self._ClassName(), "dynamic relaxation factor reaches lower bound: {}".format(self.alpha_min))
            delta_x = alpha * self.R[0]
            self.alpha_old = alpha
            return delta_x

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "init_alpha"     :  0.1,
            "init_alpha_max" :  0.45,
            "alpha_max"      :  2.0,
            "alpha_min"      : -2.0
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
