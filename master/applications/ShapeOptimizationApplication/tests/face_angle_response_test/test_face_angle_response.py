import math

# Import Kratos core and apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase
from KratosMultiphysics.ShapeOptimizationApplication.response_functions.face_angle import FaceAngleResponseFunction


class FaceAngleTest(TestCase):

    def _test(self, settings, ref_gradient, ref_value):
        model = KM.Model()
        mp = model.CreateModelPart("surface")

        response = FaceAngleResponseFunction("face_angle", settings, model)

        mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

        mp.AddNodalSolutionStepVariable(KSO.SHAPE_UPDATE)
        mp.AddNodalSolutionStepVariable(KSO.SHAPE_CHANGE)
        mp.AddNodalSolutionStepVariable(KM.NORMAL)
        mp.AddNodalSolutionStepVariable(KSO.NORMALIZED_SURFACE_NORMAL)

        mp.CreateNewNode(1, 1.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 1.0, 0.0)
        mp.CreateNewNode(3, 0.0, 0.0, 1.0)
        mp.CreateNewNode(4, 0.0, 1.0, 1.0)
        mp.CreateNewNode(5, 0.0, 0.0, 2.0)
        mp.CreateNewNode(6, 0.0, 1.0, 2.0)
        mp.CreateNewNode(7, 1.0, 0.0, 3.0)
        mp.CreateNewNode(8, 1.0, 1.0, 3.0)

        prop = mp.GetProperties()[1]

        mp.CreateNewCondition("SurfaceCondition3D4N", 1, [1, 2, 4, 3], prop)
        mp.CreateNewCondition("SurfaceCondition3D4N", 2, [3, 4, 6, 5], prop)
        mp.CreateNewCondition("SurfaceCondition3D4N", 3, [5, 6, 8, 7], prop)

        response.Initialize()
        response.InitializeSolutionStep()
        response.CalculateValue()
        response.CalculateGradient()
        value = response.GetValue()
        gradient = response.GetNodalGradient(KM.SHAPE_SENSITIVITY)

        self.assertAlmostEqual(value, ref_value)

        for node in mp.Nodes:
            node_id = node.Id
            self.assertVectorAlmostEqual(KM.Vector(gradient[node_id]), KM.Vector(ref_gradient[node_id]))

    def test_min_angle_0(self):
        settings = KM.Parameters("""{
            "response_type"         : "face_angle",
            "model_part_name"       : "surface",
            "model_import_settings" : {
                "input_type"        : "use_input_model_part"
            },
            "main_direction": [0.0, 0.0, 1.0],
            "min_angle": 0.0
        }""")

        ref_gradient = {
            1: [0.0, 0.0, 0.0],
            2: [0.0, 0.0, 0.0],
            3: [0.0, 0.0, 0.0],
            4: [0.0, 0.0, 0.0],
            5: [-0.1767768057492347, 0.0, 0.17677667329962787],
            6: [-0.1767768057492347, 0.0, 0.17677667329962787],
            7: [0.1767765847038305, 0.0, -0.17677671737548195],
            8: [0.1767765848148528, 0.0, -0.17677671737548195],
        }

        ref_value = 0.7071067811865475

        self._test(settings, ref_gradient, ref_value)

    def test_min_angle_5(self):
        settings = KM.Parameters("""{
            "response_type"         : "face_angle",
            "model_part_name"       : "surface",
            "model_import_settings" : {
                "input_type"        : "use_input_model_part"
            },
            "main_direction": [0.0, 0.0, 1.0],
            "min_angle": 5.0
        }""")

        ref_gradient = {
            1: [0.0, 0.0, 0.0],
            2: [0.0, 0.0, 0.0],
            3: [-0.054538461149566224, 0.0, 0.0],
            4: [-0.054538461149566224, 0.0, 0.0],
            5: [-0.12118357252354553, 0.0, 0.17572190201379043],
            6: [-0.12118357252354553, 0.0, 0.17572190201379043],
            7: [0.17572181394661654, 0.0, -0.17572194582665762],
            8: [0.17572181405697643, 0.0, -0.17572194582665762],
        }

        ref_value = 0.7990300873059978

        self._test(settings, ref_gradient, ref_value)

    def test_min_angle_minus_5(self):
        settings = KM.Parameters("""{
            "response_type"         : "face_angle",
            "model_part_name"       : "surface",
            "model_import_settings" : {
                "input_type"        : "use_input_model_part"
            },
            "main_direction": [0.0, 0.0, 1.0],
            "min_angle": -5.0
        }""")

        ref_gradient = {
            1: [0.0, 0.0, 0.0],
            2: [0.0, 0.0, 0.0],
            3: [0.0, 0.0, 0.0],
            4: [0.0, 0.0, 0.0],
            5: [-0.1767768057492347, 0.0, 0.17677667329962787],
            6: [-0.1767768057492347, 0.0, 0.17677667329962787],
            7: [0.1767765847038305, 0.0, -0.17677671737548195],
            8: [0.1767765848148528, 0.0, -0.17677671737548195],
        }

        ref_value = 0.6199510384388893

        self._test(settings, ref_gradient, ref_value)


if __name__ == '__main__':
    KM.KratosUnittest.main()
