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

    def test_rotating_frame_process(self):
        model, model_part = GenerateModel()

        self.print_reference_values = False

        parameters = KratosMultiphysics.Parameters("""
        {
            "rotating_frame_model_part_name": "Main.GENERIC_RotatingFrame",
            "rotating_object_model_part_name": "Main.GENERIC_RotatingObject",
            "center_of_rotation": [0.0, 0.0, 0.0],
            "axis_of_rotation": [0.0, 0.0, 1.0],
            "target_angular_velocity_radians": -10.0,
            "acceleration_time": 0.0,
            "fix_mesh_displacement": true,
            "fix_velocity": true
        }""")

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


if __name__ == "__main__":
    KratosUnittest.main()