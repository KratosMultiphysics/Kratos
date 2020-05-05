import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.ChimeraApplication as chm
from KratosMultiphysics.ChimeraApplication.rotate_region_process import ApplyRotateRegionProcess

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
                "model_part_name":"Main",
                "center_of_rotation":[0.0,0.0,0.0],
                "angular_velocity_radians": 1.5708,
                "axis_of_rotation":[1.0,1.0,1.0],
                "is_ale" : false
            }""")
        rotate_region_process = ApplyRotateRegionProcess(current_model, rotation_parameters)

        rotate_region_process.ExecuteBeforeSolutionLoop()
        rotate_region_process.ExecuteInitializeSolutionStep()
        rotate_region_process.ExecuteFinalizeSolutionStep()
        rotate_region_process.ExecuteAfterOutputStep()

        self.__CheckRotation_only_rotation(model_part)

    def test_RotateRegionProcessWithTorque(self):
        current_model = KratosMultiphysics.Model()
        model_part_name = "Main"
        model_part = current_model.CreateModelPart(model_part_name)
        self.__MakeModelPart(model_part)
        model_part.CloneTimeStep(1.0)

        rotation_parameters = KratosMultiphysics.Parameters("""{
                "model_part_name":"Main",
                "center_of_rotation":[3.0,2.5,2.5],
				"calculate_torque" : true,
				"moment_of_inertia" : 12,
				"rotational_damping" :2,
                "axis_of_rotation":[1.0,0.0,0.0],
                "is_ale" : false
            }""")
        rotate_region_process = ApplyRotateRegionProcess(current_model, rotation_parameters)

        reaction_one = 1.0
        reaction_two = 1.0
        time = 0.0
        dt = 0.2
        step = 0
        end_time = 10.0
        rotate_region_process.ExecuteInitialize()
        rotate_region_process.ExecuteBeforeSolutionLoop()

        while (time <= end_time):
            time = time + dt
            step = step + 1
            model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
            rotate_region_process.ExecuteInitializeSolutionStep()

            model_part.CloneTimeStep(time)
            model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.REACTION_X, 0, reaction_one)
            model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.REACTION_Y, 0, reaction_one)
            model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.REACTION_Z, 0, reaction_one)

            model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.REACTION_X, 0, -reaction_two)
            model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.REACTION_Y, 0, -reaction_two)
            model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.REACTION_Z, 0, -reaction_two)
            rotate_region_process.ExecuteFinalizeSolutionStep()

        rotate_region_process.ExecuteAfterOutputStep()
        rotate_region_process.ExecuteFinalize()
        self.__CheckRotation_with_torque(model_part)

    def __MakeModelPart(self, model_part):
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        model_part.AddNodalSolutionStepVariable(chm.ROTATION_MESH_DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(chm.ROTATION_MESH_VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)

        # Create nodes
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        model_part.CreateNewNode(1, 0.00000, 1.00000, 2.00000)
        model_part.CreateNewNode(2, 0.00000, 0.50000, 1.00000)
        model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DENSITY,1)
        model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DENSITY,1)

    def __CheckRotation_only_rotation(self, model_part):
        self.assertAlmostEqual(model_part.Nodes[1].X, 1.5773539423908345)
        self.assertAlmostEqual(model_part.Nodes[1].Y, -0.1547005383714617)
        self.assertAlmostEqual(model_part.Nodes[1].Z, 1.577346595980627)

        self.assertAlmostEqual(model_part.Nodes[2].X, 0.7886769711954172)
        self.assertAlmostEqual(model_part.Nodes[2].Y, -0.07735026918573085)
        self.assertAlmostEqual(model_part.Nodes[2].Z, 0.7886732979903135)

    def __CheckRotation_with_torque(self, model_part):
        self.assertAlmostEqual(model_part.Nodes[2].Z, 0.9637133210892261)
        self.assertAlmostEqual(model_part.GetValue(chm.ROTATIONAL_ANGLE), 1.8363209179941489)