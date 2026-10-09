# Kratos Imports
import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.json_output_process import JsonOutputProcess
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess


mesh_moving_is_available = kratos_utilities.CheckIfApplicationsAvailable("MeshMovingApplication")


def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), file_name)


def GenerateModel():
    model = KratosMultiphysics.Model()
    model_part = model.CreateModelPart("Main")

    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

    model_part_io = KratosMultiphysics.ModelPartIO(
        GetFilePath("test_mdpa_files/rotating_frame_process_test")
    )
    model_part_io.ReadModelPart(model_part)

    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MESH_DISPLACEMENT_X, model_part)
    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MESH_DISPLACEMENT_Y, model_part)
    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MESH_DISPLACEMENT_Z, model_part)
    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_X, model_part)
    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Y, model_part)
    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Z, model_part)

    model_part.ProcessInfo[KratosMultiphysics.TIME] = 1.0
    model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = 1.0

    return model, model_part


@KratosUnittest.skipUnless(mesh_moving_is_available, "MeshMovingApplication is not available")
class TestRotatingFrameProcess(KratosUnittest.TestCase):

    @staticmethod
    def _GetBaseParameters():
        return KratosMultiphysics.Parameters("""
        {
            "rotating_frame_model_part_name": "Main.GENERIC_RotatingFrame",
            "rotating_object_model_part_name": "Main.GENERIC_RotatingObject",
            "center_of_rotation": [0.0, 0.0, 0.0],
            "axis_of_rotation": [0.0, 0.0, 1.0],
            "fix_mesh_displacement": true,
            "fix_velocity": true
        }""")

    def test_rotating_frame_process(self):
        model, model_part = GenerateModel()

        self.print_reference_values = False

        parameters = self._GetBaseParameters()
        parameters.AddEmptyValue("angular_velocity_radians").SetDouble(-10.0)

        from KratosMultiphysics.MeshMovingApplication.rotating_frame_process import RotatingFrameProcess

        process = RotatingFrameProcess(model, parameters)
        process.ExecuteInitializeSolutionStep()

        if self.print_reference_values:
            json_output_settings = KratosMultiphysics.Parameters(r'''{
                "output_variables": ["VELOCITY", "MESH_DISPLACEMENT"],
                "output_file_name": "",
                "model_part_name": "Main",
                "time_frequency": 0.0
            }''')
            json_output_settings["output_file_name"].SetString(
                GetFilePath("test_results/rotating_frame_process_test_results.json")
            )

            json_output_process = JsonOutputProcess(model, json_output_settings)
            json_output_process.ExecuteInitialize()
            json_output_process.ExecuteBeforeSolutionLoop()
            json_output_process.ExecuteFinalizeSolutionStep()
        else:
            json_check_parameters = KratosMultiphysics.Parameters(r'''{
                "check_variables": ["VELOCITY", "MESH_DISPLACEMENT"],
                "input_file_name": "",
                "model_part_name": "Main",
                "time_frequency": 0.0,
                "tolerance": 1e-6
            }''')
            json_check_parameters["input_file_name"].SetString(
                GetFilePath("test_results/rotating_frame_process_test_results.json")
            )

            json_check_process = FromJsonCheckResultProcess(model, json_check_parameters)
            json_check_process.ExecuteInitialize()
            json_check_process.ExecuteBeforeSolutionLoop()
            json_check_process.ExecuteFinalizeSolutionStep()

    def test_rotating_frame_process_with_angular_velocity_string(self):
        numeric_model, numeric_model_part = GenerateModel()
        function_model, function_model_part = GenerateModel()

        from KratosMultiphysics.MeshMovingApplication.rotating_frame_process import RotatingFrameProcess

        numeric_parameters = self._GetBaseParameters()
        numeric_parameters.AddEmptyValue("angular_velocity_radians").SetDouble(-10.0)

        function_parameters = self._GetBaseParameters()
        function_parameters.AddEmptyValue("angular_velocity_radians").SetString("-10.0")

        numeric_process = RotatingFrameProcess(numeric_model, numeric_parameters)
        function_process = RotatingFrameProcess(function_model, function_parameters)

        numeric_process.ExecuteInitializeSolutionStep()
        function_process.ExecuteInitializeSolutionStep()

        for numeric_node in numeric_model_part.Nodes:
            function_node = function_model_part.GetNode(numeric_node.Id)
            self.assertVectorAlmostEqual(
                numeric_node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT),
                function_node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT),
                8
            )
            self.assertVectorAlmostEqual(
                numeric_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY),
                function_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY),
                8
            )

    def test_rotating_frame_process_outside_interval_is_no_op(self):
        model, model_part = GenerateModel()
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 1.0

        from KratosMultiphysics.MeshMovingApplication.rotating_frame_process import RotatingFrameProcess

        parameters = self._GetBaseParameters()
        parameters.AddValue("interval", KratosMultiphysics.Parameters('[2.0, "End"]'))
        parameters.AddEmptyValue("angular_velocity_radians").SetDouble(-10.0)

        process = RotatingFrameProcess(model, parameters)
        process.ExecuteInitializeSolutionStep()

        for node in model_part.Nodes:
            self.assertVectorAlmostEqual(
                node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT),
                [0.0, 0.0, 0.0],
                12
            )
            self.assertVectorAlmostEqual(
                node.GetSolutionStepValue(KratosMultiphysics.VELOCITY),
                [0.0, 0.0, 0.0],
                12
            )

    def test_rotating_frame_process_split_intervals_match_piecewise_function(self):
        piecewise_model, piecewise_model_part = GenerateModel()
        split_model, split_model_part = GenerateModel()

        from KratosMultiphysics.MeshMovingApplication.rotating_frame_process import RotatingFrameProcess

        piecewise_parameters = self._GetBaseParameters()
        piecewise_parameters["fix_mesh_displacement"].SetBool(False)
        piecewise_parameters["fix_velocity"].SetBool(False)
        piecewise_parameters.AddEmptyValue("angular_velocity_radians").SetString("-10.0*t if t<1.0 else -10.0")

        split_parameters_1 = self._GetBaseParameters()
        split_parameters_1["fix_mesh_displacement"].SetBool(False)
        split_parameters_1["fix_velocity"].SetBool(False)
        split_parameters_1.AddValue("interval", KratosMultiphysics.Parameters('[0.0, 1.0]'))
        split_parameters_1.AddEmptyValue("angular_velocity_radians").SetString("-10.0*t")

        split_parameters_2 = self._GetBaseParameters()
        split_parameters_2["fix_mesh_displacement"].SetBool(False)
        split_parameters_2["fix_velocity"].SetBool(False)
        split_parameters_2.AddValue("interval", KratosMultiphysics.Parameters('[1.0, "End"]'))
        split_parameters_2.AddEmptyValue("angular_velocity_radians").SetDouble(-10.0)

        piecewise_process = RotatingFrameProcess(piecewise_model, piecewise_parameters)
        split_process_1 = RotatingFrameProcess(split_model, split_parameters_1)
        split_process_2 = RotatingFrameProcess(split_model, split_parameters_2)

        # Step at t=1.0 (dt=1.0)
        piecewise_process.ExecuteInitializeSolutionStep()
        split_process_1.ExecuteInitializeSolutionStep()
        split_process_2.ExecuteInitializeSolutionStep()

        # Step at t=2.0 (dt=1.0)
        piecewise_model_part.ProcessInfo[KratosMultiphysics.TIME] = 2.0
        split_model_part.ProcessInfo[KratosMultiphysics.TIME] = 2.0
        piecewise_process.ExecuteInitializeSolutionStep()
        split_process_1.ExecuteInitializeSolutionStep()
        split_process_2.ExecuteInitializeSolutionStep()

        for piecewise_node in piecewise_model_part.Nodes:
            split_node = split_model_part.GetNode(piecewise_node.Id)
            self.assertVectorAlmostEqual(
                piecewise_node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT),
                split_node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT),
                8
            )
            self.assertVectorAlmostEqual(
                piecewise_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY),
                split_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY),
                8
            )


if __name__ == "__main__":
    KratosUnittest.main()
