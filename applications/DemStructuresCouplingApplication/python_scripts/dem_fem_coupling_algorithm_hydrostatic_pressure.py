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
        return
        '''self.initial_rigid_face_model_part = self.dem_solution.rigid_face_model_part
        self.initial_rigid_face_area = 0.0
        for elem in self.initial_rigid_face_model_part.Conditions:
            self.initial_rigid_face_area += elem.GetGeometry().Area()
        print(self.initial_rigid_face_area)'''

        self.top_fem_model_part = self.model["Structure.DISPLACEMENT_toppressure"]
        self.top_fem_model_part_area = 0.0
        for elem in self.top_fem_model_part.Conditions:
            self.top_fem_model_part_area += elem.GetGeometry().Area()
        print(self.top_fem_model_part_area)

        '''num_of_nodes = 0
        for node in self.initial_rigid_face_model_part.Nodes:
            if node.Y > 0.003:
                num_of_nodes += 1
        print(num_of_nodes)'''

    def PreviousCalculations(self):
        return

        for node in self.dem_solution.spheres_model_part.Nodes:
            if node.Y > 0.003:
                node.SetSolutionStepValue(Kratos.VELOCITY_Y, -0.3)
                node.Fix(Kratos.VELOCITY_Y)

        for node in self.dem_solution.spheres_model_part.Nodes:
            if node.Y < 0.001:
                node.SetSolutionStepValue(Kratos.VELOCITY_Y, 0.3)
                node.Fix(Kratos.VELOCITY_Y)

        reaction_forces_y = 0.0
        total_reaction_forces_y = 0.0

        for node in self.dem_solution.spheres_model_part.Nodes:
            if node.Y < 0.001:
                reaction_forces_y = node.GetSolutionStepValue(Dem.FORCE_REACTION_Y)
                total_reaction_forces_y += reaction_forces_y

        #Octogon area = 0.000513
        octogon_area = 0.000513

        sigma_DEM_in_Pa = total_reaction_forces_y / octogon_area

        print("sigma_DEM_in_Pa")
        print(sigma_DEM_in_Pa)

        '''num_of_nodes = 0
        elastic_forces_y = 0.0
        total_elastic_forces_y = 0.0

        for node in self.initial_rigid_face_model_part.Nodes:
            if node.Y > 0.003:
                elastic_forces_y = node.GetSolutionStepValue(Dem.ELASTIC_FORCES_Y)
                total_elastic_forces_y += elastic_forces_y
                if elastic_forces_y:
                    num_of_nodes += 1
        print(num_of_nodes)'''

        num_of_nodes = 0
        reaction_forces_y = 0.0
        total_reaction_forces_y = 0.0
        for node in self.top_fem_model_part.Nodes:
            reaction_forces_y = node.GetSolutionStepValue(Kratos.REACTION_Y)
            total_reaction_forces_y += reaction_forces_y

        sigma_FEM_in_Pa = -total_reaction_forces_y / self.top_fem_model_part_area
        print("sigma_FEM_in_Pa")
        print(sigma_FEM_in_Pa)

        '''self.top_pressure_model_part = self.model["Structure.DISPLACEMENT_toppressure"]
        self.bottom_pressure_model_part = self.model["Structure.DISPLACEMENT_bottompressure"]'''
        self.outer_walls_model_part = self.model["Structure.SurfacePressure3D_lateralpressure"]

        variable_utils = Kratos.VariableUtils()
        '''variable_utils.SetScalarVar(Kratos.POSITIVE_FACE_PRESSURE, sigma_DEM_in_Pa, self.top_pressure_model_part.Nodes)
        variable_utils.SetScalarVar(Kratos.POSITIVE_FACE_PRESSURE, sigma_DEM_in_Pa, self.bottom_pressure_model_part.Nodes)'''
        variable_utils.SetScalarVar(Kratos.POSITIVE_FACE_PRESSURE, sigma_DEM_in_Pa, self.outer_walls_model_part.Nodes)

if __name__ == "__main__":
    Algorithm().Run()
