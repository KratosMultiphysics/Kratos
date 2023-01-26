# This script contains an algorithm that models fluid-particle interaction.
# It combines two parts: a FEM model for the fluid and a DEM model for the particles.
# It has been conceived by adding the DEM part and the interaction on top of an original fluid-only script (see kratos/applications/FluidDynamicswApplication)
# Some parts of the original fluid script have been kept practically untouched and are clearly marked.
# Whenever a minor modification has been made on one of these parts, the corresponding line is indicated with a comment: # MOD.

import sys

# Kratos
import KratosMultiphysics as KratosMultiphysics
from KratosMultiphysics import Model, Parameters
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.IncompressibleFluidApplication
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.SwimmingDEMApplication
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.PfemSolidMechanicsApplication
import KratosMultiphysics.PfemFluidDynamicsApplication
import KratosMultiphysics.ContactMechanicsApplication
import KratosMultiphysics.ExternalSolversApplication


class Solution:
    def __init__(self, model, algorithm = None, varying_parameters = Parameters("{}")):
        self.alg = algorithm

        if self.alg == None:
            import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_PFEM_algorithm as swimming_DEM_PFEM_algorithm
            self.alg = swimming_DEM_PFEM_algorithm.Algorithm(model, varying_parameters)

    def Run(self):
        return self.alg.Run()

if __name__=="__main__":
    model = Model()
    Solution(model).Run()
