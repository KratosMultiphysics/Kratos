## @module aitken
# This module contains the class Aitken
# Author: Wei He
# Date: Feb. 20, 2017
try :
    import numpy as np
except ModuleNotFoundError:
    print(tools.bcolors.FAIL + 'Numpy is not available ! '+ tools.bcolors.ENDC)
    exit()

from copy import deepcopy
from collections import deque

## Class Aitken.
# This class contains the implementation of Aitken relaxation and helper functions.
# Reference: Ulrich Küttler et al., "Fixed-point fluid–structure interaction solvers with dynamic relaxation"
class Aitken:
    ## The constructor.
    # @param init_alpha Initial relaxation factor in the first time step.
    # @param init_alpha_max Maximum relaxation factor for the first iteration in each time step
    def __init__( self, init_alpha = 0.1, init_alpha_max = 0.45 ):
            self.R = deque( maxlen = 2 )
            self.init_alpha_max = init_alpha_max
            self.alpha_old = init_alpha
            print("Quasi-Newton method: Aitken relaxation")

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
                print( "Aitken: Doing relaxation in the first iteration with initial factor = ", alpha )
                return alpha * r
            else:
                r_diff = self.R[0] - self.R[1]
                numerator = np.inner( self.R[1], r_diff )
                denominator = np.inner( r_diff, r_diff )
                alpha = -self.alpha_old * numerator/denominator
                print( "Aitken: Doing relaxation with factor = ", alpha )
                if alpha > 20:
                   alpha = 20
                   print( "WARNING: dynamic relaxation factor reaches upper bound: 20" )
                elif alpha < -2:
                   alpha = -2
                   print( "WARNING: dynamic relaxation factor reaches lower bound: -2" )
                delta_x = alpha * self.R[0]
            self.alpha_old = alpha

            return delta_x

    ## AdvanceTimeStep()
    # Finalizes the current time step and initializes the next time step.
    def AdvanceTimeStep( self ):
            print( "" )   # Do nothing for Aitken relaxation
