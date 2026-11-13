from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import time_series_inference_process
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestTimeSeriesInferenceProcess(KratosUnittest.TestCase):
    """A scripted linear extrapolator (next = 2*last - previous) is exact for
    fields linear in time, so predictions are analytically checkable."""

    def setUp(self):
        self.checkpoint = Path("test_time_series_extrapolator.pt")

        class LinearExtrapolator(torch.nn.Module):
            def forward(self, window):
                # window: (N, 2*W), oldest first -> next = 2*last - previous
                width = window.shape[1] // 2
                return 2.0 * window[:, width:] - window[:, :width]

        torch.jit.script(LinearExtrapolator()).save(str(self.checkpoint))

        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        for i in range(4):
            self.model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def _SetLinearInTime(self, t):
        # u(x, t) = (1 + x) * t : linear in time, node-dependent slope.
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, (1.0 + node.X) * t)

    def _CreateProcess(self, history_size=2):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {
                    "checkpoint_file" : "test_time_series_extrapolator.pt",
                    "device"          : "cpu"
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "history_size"    : 2
            }
        }""")
        settings["Parameters"]["history_size"].SetInt(history_size)
        return time_series_inference_process.Factory(settings, self.model)

    def test_WarmUpThenExactPrediction(self):
        process = self._CreateProcess()

        # Step 1: history 1/2 -> warm-up, no prediction.
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        self._SetLinearInTime(1.0)
        process.ExecuteFinalizeSolutionStep()
        for node in self.model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(Kratos.TEMPERATURE), 0.0)

        # Step 2: history full -> prediction of t=3 state from t=1, t=2.
        self.model_part.ProcessInfo[Kratos.STEP] = 2
        self._SetLinearInTime(2.0)
        process.ExecuteFinalizeSolutionStep()
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE), (1.0 + node.X) * 3.0, places=12)

        # Step 3: rolling window (t=2, t=3-state) -> predicts t=4.
        self.model_part.ProcessInfo[Kratos.STEP] = 3
        self._SetLinearInTime(3.0)
        process.ExecuteFinalizeSolutionStep()
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE), (1.0 + node.X) * 4.0, places=12)

    def test_OutputNormalizationFromTheCardIsApplied(self):
        """A model trained on standardized next states predicts standardized
        next states; the card's "output_normalization" makes them physical.
        The extrapolator is affine in its window, so a missing shift or
        scale cannot hide."""
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
        mean, std = 7.0, 0.5
        model_registry.SaveModelCard(self.checkpoint, {
            "output_normalization": {"type": "mean_std", "mean": [mean], "std": [std]}})
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")

        process = self._CreateProcess()
        for step, t in ((1, 1.0), (2, 2.0)):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            self._SetLinearInTime(t)
            process.ExecuteFinalizeSolutionStep()
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE),
                std * (1.0 + node.X) * 3.0 + mean, places=12)

    def test_InputNormalizationFromTheCardIsApplied(self):
        # the whole (N, K*W) window is standardized: next = 2 l' - p' with
        # l' = (l - m)/s, p' = (p - m)/s, i.e. (2 l - p - m)/s
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
        mean, std = 3.0, 2.0
        model_registry.SaveModelCard(self.checkpoint, {
            "input_normalization": {"type": "mean_std", "mean": [mean] * 2, "std": [std] * 2}})
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")

        process = self._CreateProcess()
        for step, t in ((1, 1.0), (2, 2.0)):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            self._SetLinearInTime(t)
            process.ExecuteFinalizeSolutionStep()
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE),
                ((1.0 + node.X) * 3.0 - mean) / std, places=12)

    def test_IntervalGatingSkipsHistoryToo(self):
        process = self._CreateProcess()
        settings_interval = 2
        process.output_interval = settings_interval

        for step, t in ((1, 1.0), (2, 2.0), (3, 3.0), (4, 4.0)):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            self._SetLinearInTime(t)
            process.ExecuteFinalizeSolutionStep()
        # Sampled states: t=2 (warm-up) and t=4 -> extrapolates t=6.
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE), (1.0 + node.X) * 6.0, places=12)

    def test_TooSmallHistoryRejected(self):
        with self.assertRaisesRegex(ValueError, "history_size"):
            self._CreateProcess(history_size=1)


if __name__ == '__main__':
    KratosUnittest.main()
