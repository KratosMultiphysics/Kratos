import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.ChimeraApplication import RotateRegionProcess

class ChimeraRotateRegionTest(UnitTest.TestCase):
    def setUp(self):
        pass

    def test_RotateRegionProcess(self):
        current_model = KratosMultiphysics.Model()
        model_part_name = "Main"
        model_part = current_model.CreateModelPart(model_part_name)
        self.__MakeModelPart(model_part)
        model_part.CloneTimeStep(1.0)

        rotation_parameters = KratosMultiphysics.Parameters("""{
                "center_of_rotation":[0.0,0.0,0.0],
                "angular_velocity_radians": 1.5708,
                "axis_of_rotation":[1.0,1.0,1.0],
                "is_ale" : false
            }""")
        rotate_region_process = RotateRegionProcess(model_part, rotation_parameters)

        rotate_region_process.ExecuteBeforeSolutionLoop()
        rotate_region_process.ExecuteInitializeSolutionStep()
        rotate_region_process.ExecuteFinalizeSolutionStep()
        rotate_region_process.ExecuteAfterOutputStep()

        self.__CheckRotation(model_part)

    def __MakeModelPart(self, model_part):
        # Create nodes
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        model_part.CreateNewNode(1, 0.00000, 1.00000, 2.00000)
        model_part.CreateNewNode(2, 0.00000, 0.50000, 1.00000)

    def __CheckRotation(self, model_part):
        self.assertAlmostEqual(model_part.Nodes[1].X, 1.5773539423908345)
        self.assertAlmostEqual(model_part.Nodes[1].Y, -0.1547005383714617)
        self.assertAlmostEqual(model_part.Nodes[1].Z, 1.577346595980627)

        self.assertAlmostEqual(model_part.Nodes[2].X, 0.7886769711954172)
        self.assertAlmostEqual(model_part.Nodes[2].Y, -0.07735026918573085)
        self.assertAlmostEqual(model_part.Nodes[2].Z, 0.7886732979903135)