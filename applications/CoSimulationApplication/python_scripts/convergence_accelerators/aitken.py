from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

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
        super(AitkenConvergenceAccelerator, self).__init__(settings)
        self.R = deque( maxlen = 2 )

        self.alpha_old = self.settings["init_alpha"].GetDouble()
        self.init_alpha_max = self.settings["init_alpha_max"].GetDouble()

    ## UpdateSolution(r, x)
    # @param r residual r_k
    # @param x solution x_k
    # Computes the approximated update in each iteration.
    def UpdateSolution( self, r, x ):
        self.R.appendleft( deepcopy(r) )
        k = len( self.R ) - 1
        ## For the first iteration, do relaxation only
        if k == 0:
            alpha = min( self.alpha_old, self.init_alpha_max )
            if self.echo_level > 3:
                cs_tools.cs_print_info(self._ClassName(), ": Doing relaxation in the first iteration with initial factor = " + "{0:.1g}".format(alpha))
            return alpha * r
        else:
            r_diff = self.R[0] - self.R[1]
            numerator = np.inner( self.R[1], r_diff )
            denominator = np.inner( r_diff, r_diff )
            alpha = -self.alpha_old * numerator/denominator
            if self.echo_level > 3:
                cs_tools.cs_print_info(self._ClassName(), ": Doing relaxation with factor = " + "{0:.1g}".format(alpha))
            if alpha > 20:
                alpha = 20
                if self.echo_level > 0:
                    cs_tools.cs_print_warning(self._ClassName(), "dynamic relaxation factor reaches upper bound: 20")
            elif alpha < -2:
                alpha = -2
                if self.echo_level > 0:
                    cs_tools.cs_print_warning(self._ClassName(), "dynamic relaxation factor reaches lower bound: -2")
            delta_x = alpha * self.R[0]
        self.alpha_old = alpha

        return delta_x

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "init_alpha"     : 0.1,
            "init_alpha_max" : 0.45
        }""")
        this_defaults.AddMissingParameters(super(AitkenConvergenceAccelerator, cls)._GetDefaultSettings())
        return this_defaults
