from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import math
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestCalculateDistanceToSkin(KratosUnittest.TestCase):

    def test_naca_0012_calculate_distance_to_skin_2d(self):
        # Set the problem domain using the structured mesh generator process
        model_part = KratosMultiphysics.ModelPart("ModelPart")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EMBEDDED_VELOCITY)

        problem_domain = KratosMultiphysics.Quadrilateral2D4(
            KratosMultiphysics.Node(1, -0.25, -0.75, 0.0),
            KratosMultiphysics.Node(2, -0.25,  0.75, 0.0),
            KratosMultiphysics.Node(3,  1.25,  0.75, 0.0),
            KratosMultiphysics.Node(4,  1.25, -0.75, 0.0))
        parameters = KratosMultiphysics.Parameters("{}")
        parameters.AddEmptyValue("element_name").SetString("Element2D3N")
        parameters.AddEmptyValue("condition_name").SetString("Condition2D2N")
        parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)
        parameters.AddEmptyValue("number_of_divisions").SetInt(41)

        KratosMultiphysics.StructuredMeshGeneratorProcess(problem_domain, model_part, parameters).Execute()

        # Set aerofoil geometry
        skin_model_part = KratosMultiphysics.ModelPart("Aerofoil")
        KratosMultiphysics.ModelPartIO(GetFilePath("test_calculate_distance_to_skin_naca_0012")).ReadModelPart(skin_model_part)

        # Call the CalculateDistanceToSkinProcess()
        KratosMultiphysics.CalculateDistanceToSkinProcess2D(model_part, skin_model_part).Execute()

        # Print results (left it here for debugging)
        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(
            model_part,
            "test_naca_0012_calculate_distance_to_skin_output",
            KratosMultiphysics.Parameters("""{
                "result_file_configuration": {
                    "gidpost_flags": {
                        "GiDPostMode": "GiD_PostAscii",
                        "WriteDeformedMeshFlag": "WriteUndeformed",
                        "WriteConditionsFlag": "WriteConditions",
                        "MultiFileFlag": "SingleFile"
                    },
                    "file_label": "time",
                    "output_control_type": "step",
                    "output_frequency": 1.0,
                    "body_output": true,
                    "nodal_results": ["DISTANCE"]
                }
            }"""))
        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()

    def test_naca_0012_calculate_distance_to_skin_3d(self):
        # Set the problem domain using the structured mesh generator process
        model_part = KratosMultiphysics.ModelPart("ModelPart")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EMBEDDED_VELOCITY)

        problem_domain = KratosMultiphysics.Hexahedra3D8(
            KratosMultiphysics.Node(1, -0.25, -0.75, -0.5),
            KratosMultiphysics.Node(2,  1.25, -0.75, -0.5),
            KratosMultiphysics.Node(3,  1.25,  0.75, -0.5),
            KratosMultiphysics.Node(4, -0.25,  0.75, -0.5),
            KratosMultiphysics.Node(5, -0.25, -0.75,  0.5),
            KratosMultiphysics.Node(6,  1.25, -0.75,  0.5),
            KratosMultiphysics.Node(7,  1.25,  0.75,  0.5),
            KratosMultiphysics.Node(8, -0.25,  0.75,  0.5))
        parameters = KratosMultiphysics.Parameters("{}")
        parameters.AddEmptyValue("element_name").SetString("Element3D4N")
        parameters.AddEmptyValue("condition_name").SetString("Condition3D3N")
        parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)
        parameters.AddEmptyValue("number_of_divisions").SetInt(200)

        KratosMultiphysics.StructuredMeshGeneratorProcess(problem_domain, model_part, parameters).Execute()

        # Set aerofoil geometry
        skin_model_part = KratosMultiphysics.ModelPart("Aerofoil")
        KratosMultiphysics.ModelPartIO(GetFilePath("test_calculate_distance_to_skin_naca_0012_3d")).ReadModelPart(skin_model_part)

        # Call the CalculateDistanceToSkinProcess()
        KratosMultiphysics.CalculateDistanceToSkinProcess3D(model_part, skin_model_part).Execute()

        # Print results (left it here for debugging)
        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(
            model_part,
            "test_naca_0012_calculate_distance_to_skin_3d_output",
            KratosMultiphysics.Parameters("""{
                "result_file_configuration": {
                    "gidpost_flags": {
                        "GiDPostMode": "GiD_PostAscii",
                        "WriteDeformedMeshFlag": "WriteUndeformed",
                        "WriteConditionsFlag": "WriteConditions",
                        "MultiFileFlag": "SingleFile"
                    },
                    "file_label": "time",
                    "output_control_type": "step",
                    "output_frequency": 1.0,
                    "body_output": true,
                    "nodal_results": ["DISTANCE"]
                }
            }"""))
        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()
        
if __name__ == '__main__':
    domain_size = 2
    test = TestCalculateDistanceToSkin()
    if (domain_size == 2):
        test.test_naca_0012_calculate_distance_to_skin_2d()
    else:
        test.test_naca_0012_calculate_distance_to_skin_3d()
