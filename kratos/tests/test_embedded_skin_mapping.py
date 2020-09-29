from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestEmbeddedSkinMapping(KratosUnittest.TestCase):

    def test_embedded_skin_mapping(self):
        model = KratosMultiphysics.Model()

        # Set the problem domain using the structured mesh generator process
        model_part = model.CreateModelPart("ModelPart")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        problem_domain = KratosMultiphysics.Quadrilateral2D4(
            KratosMultiphysics.Node(1, -0.25, -0.75, 0.0),
            KratosMultiphysics.Node(2, -0.25,  0.75, 0.0),
            KratosMultiphysics.Node(3,  1.25,  0.75, 0.0),
            KratosMultiphysics.Node(4,  1.25, -0.75, 0.0))
        parameters = KratosMultiphysics.Parameters("{}")
        parameters.AddEmptyValue("number_of_divisions").SetInt(100)
        parameters.AddEmptyValue("element_name").SetString("Element2D3N")
        parameters.AddEmptyValue("condition_name").SetString("LineCondition2D2N")
        parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)

        KratosMultiphysics.StructuredMeshGeneratorProcess(problem_domain, model_part, parameters).Execute()

        # Set a random field to be mapped from the model part to the skin model part
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, node.X + node.Y)

        # Set aerofoil geometry
        skin_model_part = model.CreateModelPart("SkinModelPart")
        skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_calculate_distance_to_skin_naca_0012")).ReadModelPart(skin_model_part)

        # Set a random field to be mapped from the skin model part to the (embedded) model part
        # Call the CalculateDistanceToSkinProcess()
        distance_to_skin_process = KratosMultiphysics.CalculateDistanceToSkinProcess2D(model_part, skin_model_part)
        distance_to_skin_process.Execute()

        # Set the skin intersections model part
        generated_skin_model_part = model.CreateModelPart("GeneratedSkinModelPart")
        generated_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        generated_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        embedded_skin_utility = KratosMultiphysics.EmbeddedSkinUtility2D(model_part, generated_skin_model_part, "continuous")
        embedded_skin_utility.GenerateSkin()

        # Interpolate the model part random field to the generated skin nodes
        embedded_skin_utility.InterpolateMeshVariableToSkin(KratosMultiphysics.PRESSURE, KratosMultiphysics.PRESSURE)

        # Map the pressure value from the generated intersections skins to the original skin
        # This intends emulate the step in where the loads are transferred in any coupled problem
        map_parameters = KratosMultiphysics.Parameters("""{
            "echo_level"                       : 0,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2
        }""")
        mortar_mapping_double = KratosMultiphysics.SimpleMortarMapperProcess2D2NDouble(
            generated_skin_model_part,
            skin_model_part,
            KratosMultiphysics.PRESSURE,
            KratosMultiphysics.PRESSURE,
            map_parameters)
        mortar_mapping_double.Execute()

        # Check results
        self.assertAlmostEqual(skin_model_part.GetNode(100).GetSolutionStepValue(KratosMultiphysics.PRESSURE),  0.4507103394114805)
        self.assertAlmostEqual(skin_model_part.GetNode(150).GetSolutionStepValue(KratosMultiphysics.PRESSURE),  0.1939745284224017)
        self.assertAlmostEqual(skin_model_part.GetNode(200).GetSolutionStepValue(KratosMultiphysics.PRESSURE), -0.006698755223763724)
        self.assertAlmostEqual(generated_skin_model_part.GetNode(10202).GetSolutionStepValue(KratosMultiphysics.PRESSURE), -0.002595311391968201)
        self.assertAlmostEqual(generated_skin_model_part.GetNode(10300).GetSolutionStepValue(KratosMultiphysics.PRESSURE),  0.10611968493398256)
        self.assertAlmostEqual(generated_skin_model_part.GetNode(10400).GetSolutionStepValue(KratosMultiphysics.PRESSURE),  0.3799232989903604)

if __name__ == '__main__':
    test = TestEmbeddedSkinMapping()
    test.test_embedded_skin_mapping()
