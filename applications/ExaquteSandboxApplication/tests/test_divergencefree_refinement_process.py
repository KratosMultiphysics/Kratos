import KratosMultiphysics
from KratosMultiphysics.ExaquteSandboxApplication.divergencefree_refinement_process import DivergenceFreeRefinementProcess
import KratosMultiphysics.kratos_utilities
meshing_is_available = KratosMultiphysics.kratos_utilities.CheckIfApplicationsAvailable("MeshingApplication")

import KratosMultiphysics.KratosUnittest as UnitTest


class TimeAveragingProcessTests(UnitTest.TestCase):

    def testDivergenceFreeRefinementProcess(self):
        if not meshing_is_available:
            self.skipTest("Missing required application: MeshingApplication")

        self.__CreateModel()

        settings = KratosMultiphysics.Parameters(r"""
            {
                "model_part_name": "test",
                "integration_end_point_control_value": 5.0,
                "integration_start_point_control_value": 1.0,
                "average_settings": {
                    "model_part_name": "test",
                    "variables_list": ["DIVERGENCE"],
                    "averaged_variables_list": ["AVERAGED_DIVERGENCE"],
                    "time_averaging_container": "ElementalNonHistorical",
                    "time_averaging_method": "RootMeanSquare",
                    "integration_start_point_control_variable_name": "TIME",
                    "integration_start_point_control_value": 2.0
                },
                "refinement_settings": {
                    "refinement_strategy": "global_tolerance_strategy",
                    "reference_variable_name": "AVERAGED_DIVERGENCE",
                    "minimal_size": 0.001,
                    "maximal_size": 10.0,
                    "mean_distribution_strategy": {
                        "target_refinement_coefficient": 0.9,
                        "refinement_bound": 2.0,
                        "reference_norm_name": "VELOCITY_H1_SEMINORM"
                    },
                    "maximum_strategy": {
                        "target_refinement_coefficient": 0.1,
                        "refinement_coefficient": 2.0
                    },
                    "global_tolerance_strategy": {
                        "global_tolerance": 0.1
                    }
                }
            }""")

        self.process = DivergenceFreeRefinementProcess(self.model,settings)

        velocity_vector = []
        for _ in self.model.GetModelPart("test").Nodes:
            velocity_vector.append([0.0, 0.0, 0.0])

        self.process.Check()
        self.process.ExecuteInitialize()

        current_time = 0.0
        for _ in range(1, 10):
            current_time += 0.5
            self.model.GetModelPart("test").CloneTimeStep(current_time)
            self.model.GetModelPart("test").ProcessInfo[KratosMultiphysics.STEP] += 1

            self.process.ExecuteInitializeSolutionStep()

            for node_index, node in enumerate(self.model.GetModelPart("test").Nodes):
                velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                temp_velocity = KratosMultiphysics.Array3()
                temp_velocity[0] = velocity[0] + current_time
                temp_velocity[1] = velocity[1] + current_time
                temp_velocity[2] = velocity[2] + current_time
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, temp_velocity)

                if (self.model.GetModelPart("test").ProcessInfo[KratosMultiphysics.TIME] >= 2.0):
                    velocity_vector[node_index][0] += temp_velocity[0] * 0.5
                    velocity_vector[node_index][1] += temp_velocity[1] * 0.5
                    velocity_vector[node_index][2] += temp_velocity[2] * 0.5

            self.process.ExecuteFinalizeSolutionStep()

        self.process.ExecuteFinalize()

        self.assertEqual(self.model.GetModelPart("test").NumberOfElements(),66317)
        self.assertEqual(self.model.GetModelPart("test").NumberOfNodes(),33514)

    def __CreateModel(self):
        self.model = KratosMultiphysics.Model()
        self.model.CreateModelPart("test")
        self.model.GetModelPart("test").AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.model.GetModelPart("test").ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2

        self.model.GetModelPart("test").CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model.GetModelPart("test").CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model.GetModelPart("test").CreateNewNode(3, 2.0, 0.0, 0.0)
        self.model.GetModelPart("test").CreateNewNode(4, 2.0, 1.0, 0.0)
        self.model.GetModelPart("test").CreateNewNode(5, 1.5, 1.0, 0.0)
        self.model.GetModelPart("test").CreateNewNode(6, 1.0, 1.0, 0.0)

        self.model.GetModelPart("test").AddProperties(KratosMultiphysics.Properties(1))
        self.model.GetModelPart("test").CreateNewElement("Element2D3N", 1, [1,2,6], self.model.GetModelPart("test").GetProperties()[1])
        self.model.GetModelPart("test").CreateNewElement("Element2D3N", 2, [2,5,6], self.model.GetModelPart("test").GetProperties()[1])
        self.model.GetModelPart("test").CreateNewElement("Element2D3N", 3, [2,3,5], self.model.GetModelPart("test").GetProperties()[1])
        self.model.GetModelPart("test").CreateNewElement("Element2D3N", 4, [3,4,5], self.model.GetModelPart("test").GetProperties()[1])

        coefficient = 1.0
        for node in self.model.GetModelPart("test").Nodes:
            vector = KratosMultiphysics.Vector(3)
            vector[0] = coefficient*1
            vector[1] = coefficient*2
            vector[2] = coefficient*3
            if (node.SolutionStepsDataHas(KratosMultiphysics.VELOCITY)):
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, vector)
            coefficient += coefficient

if __name__ == '__main__':
    UnitTest.main()
