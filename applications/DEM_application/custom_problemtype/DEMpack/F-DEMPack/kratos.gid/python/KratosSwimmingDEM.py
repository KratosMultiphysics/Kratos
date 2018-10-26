# This script contains an algorithm that models fluid-particle interaction.
# It combines two parts: a FEM model for the fluid and a DEM model for the particles.
# It has been conceived by adding the DEM part and the interaction on top of an original fluid-only script (see kratos/applications/FluidDynamicswApplication)
# Some parts of the original fluid script have been kept practically untouched and are clearly marked.
# Whenever a minor modification has been made on one of these parts, the corresponding line is indicated with a comment: # MOD.

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import sys

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication   import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

class Solution:
    def __enter__ (self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        pass

    def __init__(self, model, algorithm = None, varying_parameters = Parameters("{}")):

        if algorithm == None:
            import swimming_DEM_algorithm
            self.alg = swimming_DEM_algorithm.Algorithm(model, varying_parameters)
        else:
            self.alg = algorithm.Algorithm(model, varying_parameters)

    def Run(self):
        return self.alg.Run()

if __name__=="__main__":
    model = Model()
    Solution(model).Run()
