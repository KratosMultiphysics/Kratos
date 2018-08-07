# This script contains an algorithm that models fluid-particle interaction.
# It combines two parts: a FEM model for the fluid and a DEM model for the particles.
# It has been conceived by adding the DEM part and the interaction on top of an original fluid-only script (see kratos/applications/FluidDynamicswApplication)
# Some parts of the original fluid script have been kept practically untouched and are clearly marked.
# Whenever a minor modification has been made on one of these parts, the corresponding line is indicated with a comment: # MOD.

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import sys

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.DelaunayMeshingApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *
from KratosMultiphysics.PfemFluidDynamicsApplication import *
from KratosMultiphysics.ContactMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *


class Solution:
    def __init__(self, algorithm = None, varying_parameters = Parameters("{}")):
        self.alg = algorithm

        if self.alg == None:
            import swimming_DEM_PFEM_algorithm
            self.alg = swimming_DEM_PFEM_algorithm.Algorithm(varying_parameters)

    def Run(self):
        return self.alg.Run()

if __name__=="__main__":
    Solution().Run()
