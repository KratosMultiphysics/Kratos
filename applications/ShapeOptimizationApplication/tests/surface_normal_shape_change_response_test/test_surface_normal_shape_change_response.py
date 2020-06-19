import math

# Import Kratos core and apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase
from KratosMultiphysics.ShapeOptimizationApplication.response_functions.surface_normal_shape_change import SurfaceNormalShapeChange


class SurfaceNormalShapeChangeTest(TestCase):

    def test_simple_model_part(self):
        model = KM.Model()
        mp = model.CreateModelPart("surface")
        mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

        mp.AddNodalSolutionStepVariable(KSO.SHAPE_UPDATE)
        mp.AddNodalSolutionStepVariable(KSO.SHAPE_CHANGE)
        mp.AddNodalSolutionStepVariable(KM.NORMAL)
        mp.AddNodalSolutionStepVariable(KSO.NORMALIZED_SURFACE_NORMAL)

        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 1.0, 1.0, 0.0)
        mp.CreateNewNode(4, 0.0, 1.0, 0.0)
        mp.CreateNewNode(5, 1.0, 1.0, -1.0)
        mp.CreateNewNode(6, 1.0, 0.0, -1.0)

        prop = mp.GetProperties()[1]

        mp.CreateNewCondition("SurfaceCondition3D4N", 1, [1, 2, 3, 4], prop)
        mp.CreateNewCondition("SurfaceCondition3D4N", 2, [2, 6, 5, 3], prop)

        settings = KM.Parameters("""{
            "response_type"         : "surface_normal_shape_change",
            "model_part_name"       : "surface",
            "model_import_settings" : {
                "input_type"        : "use_input_model_part"
            }
        }""")

        response = SurfaceNormalShapeChange("surface_normal", settings, model)

        response.Initialize()
        response.InitializeSolutionStep()
        response.CalculateValue()
        response.CalculateGradient()
        value = response.GetValue()
        gradient = response.GetNodalGradient(KM.SHAPE_SENSITIVITY)

        ref_gradient = {
            1: [0.0, 0.0, 1.0],
            2: [math.sqrt(2.0)/2.0, 0.0, math.sqrt(2.0)/2.0],
            3: [math.sqrt(2.0)/2.0, 0.0, math.sqrt(2.0)/2.0],
            4: [0.0, 0.0, 1.0],
            5: [1.0, 0.0, 0.0],
            6: [1.0, 0.0, 0.0],
        }
        self.assertAlmostEqual(value, 0.0)
        for node in mp.Nodes:
            node_id = node.Id
            self.assertVectorAlmostEqual(KM.Vector(gradient[node_id]), KM.Vector(ref_gradient[node_id]))

        # update nodes and test again
        for node in mp.Nodes:
            node.SetSolutionStepValue(KSO.SHAPE_UPDATE, [1, 0, 0])
        KSO.MeshControllerUtilities(mp).UpdateMeshAccordingInputVariable(KSO.SHAPE_UPDATE)

        response.InitializeSolutionStep()
        response.CalculateValue()
        value = response.GetValue()
        self.assertAlmostEqual(value, 2.0 + math.sqrt(2.0))


if __name__ == '__main__':
    KM.KratosUnittest.main()
