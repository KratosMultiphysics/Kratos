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

    def PreviousCalculations(self):
        self.ComputeLoadingPath()

    def ComputeLoadingPath(self):

        self.top_pressure_model_part = self.model["Structure.DISPLACEMENT_toppressure"]
        self.bottom_pressure_model_part = self.model["Structure.DISPLACEMENT_bottompressure"]
        self.outer_walls_model_part = self.model["Structure.SurfacePressure3D_lateralpressure"]

        #Test A
        time_0 = 0.000001
        time_1 = 0.000007
        #Test B
        time_0 = 0.000001
        time_1 = 0.000002

        if self.dem_solution.time < time_0:
            sigmaZ_in_Pa = 10000 * 500 * self.dem_solution.time / 0.000145
        else:
            sigmaZ_in_Pa = 500 / 0.000145

        if self.dem_solution.time < time_0:
            sigmaY_in_Pa = 10000 * 500 * self.dem_solution.time / 0.000145
        elif self.dem_solution.time > time_0 and self.dem_solution.time < time_1:
            sigmaY_in_Pa = 500 / 0.000145 + (4000 - 500) * 10000 / (0.000145 * (time_1 - time_0)) * (self.dem_solution.time - time_0)
        else:
            sigmaY_in_Pa = 10000 * 3500 / 0.000145

        if self.dem_solution.time < time_0:
            sigmaX_in_Pa = 10000 * 500 * self.dem_solution.time / 0.000145
        else:
            sigmaX_in_Pa = 500 / 0.000145 + (4000 - 500) * 10000 / (0.000145 * (time_1 - time_0)) * (self.dem_solution.time - time_0)

        variable_utils = Kratos.VariableUtils()
        variable_utils.SetScalarVar(Kratos.POSITIVE_FACE_PRESSURE, sigmaX_in_Pa, self.top_pressure_model_part.Nodes)
        variable_utils.SetScalarVar(Kratos.POSITIVE_FACE_PRESSURE, sigmaY_in_Pa, self.bottom_pressure_model_part.Nodes)
        variable_utils.SetScalarVar(Kratos.POSITIVE_FACE_PRESSURE, sigmaZ_in_Pa, self.outer_walls_model_part.Nodes)

if __name__ == "__main__":
    Algorithm().Run()
