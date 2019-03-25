from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as Dem

import dem_fem_coupling_algorithm as CouplingAlgorithm
BaseAlgorithm = CouplingAlgorithm.Algorithm

import math

class Algorithm(BaseAlgorithm):

    def __init__(self):
        BaseAlgorithm.__init__(self)

    def PreliminaryComputations(self):

        self.area = 1.0 # Compute properly in general by printing the conditions of the corresponding surface in the mdpa
        self.bottom_fem_model_part = self.model["Structure.VELOCITY_bottomwall"]
        self.actual_sigmas = open("actual_sigmas.txt", 'w')

    def PreviousCalculations(self):

        total_reaction_forces_y = 0.0
        for node in self.bottom_fem_model_part.Nodes:
            total_reaction_forces_y += node.GetSolutionStepValue(Kratos.REACTION_Z)

        sigma_FEM_in_Pa = total_reaction_forces_y / self.area
        print("sigma_FEM_in_Pa")
        print(sigma_FEM_in_Pa)
        self.actual_sigmas.write(str("%.6g"%self.dem_solution.time).rjust(13) + " " + str("%.8g"%sigma_FEM_in_Pa).rjust(15) + '\n')
        self.actual_sigmas.flush()

    def ClosingOperations(self):
        self.actual_sigmas.close()


if __name__ == "__main__":
    Algorithm().Run()
