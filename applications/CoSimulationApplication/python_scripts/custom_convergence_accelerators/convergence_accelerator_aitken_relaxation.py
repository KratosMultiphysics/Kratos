## @module aitken
# This module contains the class Aitken
# Author: Wei He
# Updated : Aditya Ghantasala
# Date: Feb. 20, 2017
try :
    import numpy as np
except ModuleNotFoundError:
    print(tools.bcolors.FAIL + 'Numpy is not available ! '+ tools.bcolors.ENDC)
    exit()

from co_simulation_base_convergence_accelerator import CoSimulationBaseConvergenceAccelerator
import co_simulation_tools as tools
from copy import deepcopy
from collections import deque

def Create(settings):
    accelerator = AitkenAccelerator(settings)
    return accelerator


## Class Aitken.
# This class contains the implementation of Aitken relaxation and helper functions.
# Reference: Ulrich Küttler et al., "Fixed-point fluid–structure interaction solvers with dynamic relaxation"
class AitkenAccelerator(CoSimulationBaseConvergenceAccelerator):
    ## The constructor.
    # @param init_alpha Initial relaxation factor in the first time step.
    # @param init_alpha_max Maximum relaxation factor for the first iteration in each time step
    def __init__( self, settings ):
        super(AitkenAccelerator, self).__init__(settings)
        default_settings = {}
        default_settings["data_list"] = list    #MANDATORY
        default_settings["settings"] = dict    #MANDATORY
        self.settings = tools.ValidateAndAssignInputParameters(default_settings, self.settings, False)
        self.R = deque( maxlen = 2 )
        self.init_alpha_max = 0.4
        self.alpha_old = self.init_alpha_max

    ## ComputeUpdate(r, x)
    # @param r residual r_k
    # @param x solution x_k
    # Computes the approximated update in each iteration.
    def ComputeUpdate( self, r, x ):
        self.R.appendleft( deepcopy(r) )
        k = len( self.R ) - 1
        ## For the first iteration, do relaxation only
        if k == 0:
            alpha = min( self.alpha_old, self.init_alpha_max )
            print( tools.bcolors.BLUE + "Aitken: Doing relaxation in the first iteration with initial factor = " + tools.bcolors.ENDC, alpha)
            return alpha * r
        else:
            r_diff = self.R[0] - self.R[1]
            numerator = np.inner( self.R[1], r_diff )
            denominator = np.inner( r_diff, r_diff )
            alpha = -self.alpha_old * numerator/denominator
            print( tools.bcolors.BLUE + "Aitken: Doing relaxation with factor = " + tools.bcolors.ENDC, alpha )
            if alpha > 2:
                alpha = 2
                print(tools.bcolors.WARNING + "WARNING: dynamic relaxation factor reaches upper bound: 2" + tools.bcolors.ENDC)
            elif alpha < -2:
                alpha = -2
                print(tools.bcolors.WARNING + "WARNING: dynamic relaxation factor reaches lower bound: -2" + tools.bcolors.ENDC)
            delta_x = alpha * self.R[0]
        self.alpha_old = alpha

        return delta_x

    ## AdvanceTimeStep()
    # Finalizes the current time step and initializes the next time step.
    def AdvanceTimeStep( self ):
        print( "" )   # Do nothing for Aitken relaxation
