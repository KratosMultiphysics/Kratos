import os
import tempfile

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.checkpoint_process import Factory


class TestCheckpointProcess(KratosUnittest.TestCase):

    def test_checkpoint_save_and_load(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            model = Kratos.Model()
            model_part = model.CreateModelPart("Main")
            model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
            node = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.TEMPERATURE, 42.0)
            model_part.ProcessInfo[Kratos.TIME] = 1.5
            model_part.ProcessInfo[Kratos.STEP] = 3

            file_name = os.path.join(temporary_directory, "checkpoint_<time>_<step>")

            def create_process(checkpoint_state):
                settings = Kratos.Parameters("""
                {
                    "Parameters": {
                        "model_part_name": "Main",
                        "checkpoint_state": "",
                        "file_name": ""
                    }
                }
                """)
                settings["Parameters"]["checkpoint_state"].SetString(checkpoint_state)
                settings["Parameters"]["file_name"].SetString(file_name)
                return Factory(settings, model)

            create_process("save").ExecuteFinalizeSolutionStep()
            node.SetSolutionStepValue(Kratos.TEMPERATURE, -1.0)
            create_process("load").ExecuteInitializeSolutionStep()

            self.assertEqual(node.GetSolutionStepValue(Kratos.TEMPERATURE), 42.0)


if __name__ == "__main__":
    KratosUnittest.main()
