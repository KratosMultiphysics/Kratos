from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import hybrid_initialization_process
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestHybridInitializationProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_hybrid_init_model.pt")

        class CallCountingTripler(torch.nn.Module):
            # affine on purpose: a dropped offset is invisible to a model
            # with f(0) = 0
            def forward(self, x):
                return 3.0 * x + 1.0

        torch.jit.script(CallCountingTripler()).save(str(self.checkpoint))

        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        for i in range(3):
            node = self.model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.VELOCITY, [float(i), 1.0, 0.0])

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def _CreateProcess(self):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {
                    "checkpoint_file" : "test_hybrid_init_model.pt",
                    "device"          : "cpu"
                },
                "input_fields"    : [ { "variable_name" : "VELOCITY", "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "DISPLACEMENT", "data_location" : "node_historical" } ]
            }
        }""")
        return hybrid_initialization_process.Factory(settings, self.model)

    def test_WarmStartBeforeSolutionLoop(self):
        process = self._CreateProcess()
        process.ExecuteBeforeSolutionLoop()
        for node in self.model_part.Nodes:
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            displacement = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            self.assertVectorAlmostEqual(displacement, [3.0 * v + 1.0 for v in velocity])

    def test_RunsExactlyOnce(self):
        process = self._CreateProcess()
        process.ExecuteBeforeSolutionLoop()
        # Change the input; a second call must NOT recompute.
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY, [100.0, 100.0, 100.0])
        process.ExecuteBeforeSolutionLoop()
        process.ExecuteInitializeSolutionStep()
        process.ExecuteFinalizeSolutionStep()
        for i, node in enumerate(self.model_part.Nodes):
            self.assertVectorAlmostEqual(
                node.GetSolutionStepValue(Kratos.DISPLACEMENT),
                [3.0 * i + 1.0, 4.0, 1.0])

    def test_ExecutionPointSettingRejected(self):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : { "checkpoint_file" : "test_hybrid_init_model.pt" },
                "input_fields"    : [ { "variable_name" : "VELOCITY", "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "DISPLACEMENT", "data_location" : "node_historical" } ],
                "execution_point" : "finalize_solution_step"
            }
        }""")
        with self.assertRaisesRegex(ValueError, "not a valid setting"):
            hybrid_initialization_process.Factory(settings, self.model)


if __name__ == '__main__':
    KratosUnittest.main()
