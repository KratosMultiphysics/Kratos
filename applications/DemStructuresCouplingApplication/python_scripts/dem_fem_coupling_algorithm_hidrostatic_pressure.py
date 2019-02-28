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
        self.initial_rigid_face_model_part = self.dem_solution.rigid_face_model_part
        self.initial_rigid_face_area = 0.0
        for elem in self.initial_rigid_face_model_part.Conditions:
            self.initial_rigid_face_area += elem.GetGeometry().Area()

    def PreviousCalculations(self):

        total_elastic_forces_y = 0.0
        for node in self.initial_rigid_face_model_part.Nodes:
            if node.Y > 0.003:
                total_elastic_forces_y += node.GetSolutionStepValue(Dem.ELASTIC_FORCES_Y)
        sigma_DEM_in_Pa = total_elastic_forces_y / (0.5 * self.initial_rigid_face_area)

        self.top_pressure_model_part = self.model["Structure.SurfacePressure3D_toppressure"]
        self.bottom_pressure_model_part = self.model["Structure.SurfacePressure3D_bottompressure"]
        self.outer_walls_model_part = self.model["Structure.SurfacePressure3D_lateralpressure"]

        variable_utils = Kratos.VariableUtils()
        variable_utils.SetScalarVar(Kratos.POSITIVE_FACE_PRESSURE, sigma_DEM_in_Pa, self.top_pressure_model_part.Nodes)
        variable_utils.SetScalarVar(Kratos.POSITIVE_FACE_PRESSURE, sigma_DEM_in_Pa, self.bottom_pressure_model_part.Nodes)
        variable_utils.SetScalarVar(Kratos.POSITIVE_FACE_PRESSURE, sigma_DEM_in_Pa, self.outer_walls_model_part.Nodes)

if __name__ == "__main__":
    Algorithm().Run()
