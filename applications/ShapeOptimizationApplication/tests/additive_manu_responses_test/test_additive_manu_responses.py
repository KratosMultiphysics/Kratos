import math

# Import Kratos core and apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase
from KratosMultiphysics.ShapeOptimizationApplication.response_functions.am_boundary_slope_angle import AMBoundarySlopeAngleResponseFunction
from KratosMultiphysics.ShapeOptimizationApplication.response_functions.am_boundary_roughness_angle import AMBoundaryRoughnessAngleResponseFunction


class AMAnglesTest(TestCase):

    def _test_slope(self, settings, ref_gradient, ref_value):
        model = KM.Model()
        mp = model.CreateModelPart("surface")

        response = AMBoundarySlopeAngleResponseFunction("AM_boundary_slope_angle", settings, model)

        mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

        mp.AddNodalSolutionStepVariable(KSO.SHAPE_UPDATE)
        mp.AddNodalSolutionStepVariable(KSO.SHAPE_CHANGE)
        mp.AddNodalSolutionStepVariable(KM.NORMAL)
        mp.AddNodalSolutionStepVariable(KSO.NORMALIZED_SURFACE_NORMAL)

        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 1.0, 0.0)
        mp.CreateNewNode(3, 0.0, 0.0, 1.0)
        mp.CreateNewNode(4, 1.0, 1.0, 1.0)


        prop = mp.GetProperties()[1]

        mp.CreateNewCondition("SurfaceCondition3D4N", 1, [1, 2, 4, 3], prop)

        response.Initialize()
        response.InitializeSolutionStep()
        response.CalculateValue()
        response.CalculateGradient()
        value = response.GetValue()
        gradient = response.GetNodalGradient(KM.SHAPE_SENSITIVITY)

        self.assertAlmostEqual(value, ref_value)

        for node in mp.Nodes:
            node_id = node.Id
            self.assertVectorAlmostEqual(KM.Vector(gradient[node_id]), KM.Vector(ref_gradient[node_id]),5)

    def _test_roughness(self, settings, ref_gradient, ref_value):
        model = KM.Model()
        mp = model.CreateModelPart("surface")

        response = AMBoundaryRoughnessAngleResponseFunction("AM_boundary_roughness_angle", settings, model)

        mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

        mp.AddNodalSolutionStepVariable(KSO.SHAPE_UPDATE)
        mp.AddNodalSolutionStepVariable(KSO.SHAPE_CHANGE)
        mp.AddNodalSolutionStepVariable(KM.NORMAL)
        mp.AddNodalSolutionStepVariable(KSO.NORMALIZED_SURFACE_NORMAL)

        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, -1.0, 1.0, 0.0)
        mp.CreateNewNode(3, 0.0, 0.0, 1.0)
        mp.CreateNewNode(4, -1.0, 1.0, 1.0)


        prop = mp.GetProperties()[1]

        mp.CreateNewCondition("SurfaceCondition3D4N", 1, [1, 2, 4, 3], prop)

        response.Initialize()
        response.InitializeSolutionStep()
        response.CalculateValue()
        response.CalculateGradient()
        value = response.GetValue()
        gradient = response.GetNodalGradient(KM.SHAPE_SENSITIVITY)

        self.assertAlmostEqual(value, ref_value)

        for node in mp.Nodes:
            node_id = node.Id            
            self.assertVectorAlmostEqual(KM.Vector(gradient[node_id]), KM.Vector(ref_gradient[node_id]),5)            

    def test_slope_angle_45(self):
        settings = KM.Parameters("""{
            "response_type"         : "AM_boundary_slope_angle",
            "model_part_name"       : "surface",
            "model_import_settings" : {
                "input_type"        : "use_input_model_part"
            },
            "main_direction": [0.0, 1.0, 0.0],
            "min_angle": 45.0
        }""")

        ref_gradient = {
            1: [-0.883884,0.883883,0],
            2: [0.883883,-0.883884,0],
            3: [-0.883884,0.8838823,0],
            4: [0.883883,-0.883884,0]
        }

        ref_value = 0.5

        self._test_slope(settings, ref_gradient, ref_value)

    def test_roughness_angle_135(self):
        settings = KM.Parameters("""{
            "response_type"         : "AM_boundary_roughness_angle",
            "model_part_name"       : "surface",
            "model_import_settings" : {
                "input_type"        : "use_input_model_part"
            },
            "main_direction": [0.0, 1.0, 0.0],
            "max_angle": 135.0
        }""")

        ref_gradient = {
            1: [0.883883,0.883883,0],
            2: [-0.883884,-0.883884,0],
            3: [0.883883,0.883883,0.0],
            4: [-0.883884,-0.883884,0.0]
        }

        ref_value = 0.5

        self._test_roughness(settings, ref_gradient, ref_value)        

if __name__ == '__main__':
    KM.KratosUnittest.main()
