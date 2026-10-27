import json
import math
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.processes import validation_metrics_process
try:
    import physicsnemo.metrics.general.mse
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestValidationMetricsProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.output_file = Path("test_validation_metrics.json")
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        # predicted = reference + 1 on every node -> mse = 1, rmse = 1, max = 1
        for i in range(4):
            node = self.model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.TEMPERATURE, float(i))       # reference
            node.SetSolutionStepValue(Kratos.PRESSURE, float(i) + 1.0)    # predicted

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.output_file))

    def _CreateProcess(self, output_interval=1, metrics='["mse", "rmse", "max_abs_error", "wasserstein"]'):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"     : "Main",
                "list_of_comparisons" : [
                    {
                        "predicted_variable" : "PRESSURE",
                        "predicted_location" : "node_historical",
                        "reference_variable" : "TEMPERATURE",
                        "reference_location" : "node_historical",
                        "metrics"            : %s
                    }
                ],
                "output_interval"     : %d,
                "output_file"         : "test_validation_metrics.json"
            }
        }""" % (metrics, output_interval))
        return validation_metrics_process.Factory(settings, self.model)

    def test_KnownMetricValues(self):
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        self.model_part.ProcessInfo[Kratos.TIME] = 0.5
        process.ExecuteFinalizeSolutionStep()

        self.assertEqual(len(process.history), 1)
        values = process.history[0]["PRESSURE_vs_TEMPERATURE"]
        self.assertAlmostEqual(values["mse"], 1.0, places=12)
        self.assertAlmostEqual(values["rmse"], 1.0, places=12)
        self.assertAlmostEqual(values["max_abs_error"], 1.0, places=12)
        self.assertGreaterEqual(values["wasserstein"], 0.0)

    def test_IntervalGating(self):
        process = self._CreateProcess(output_interval=2)
        for step in (1, 2, 3, 4):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
        self.assertEqual([r["STEP"] for r in process.history], [2, 4])

    def test_JsonReport(self):
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteFinalize()

        with open(self.output_file) as f:
            report = json.load(f)
        self.assertEqual(len(report), 1)
        self.assertAlmostEqual(report[0]["PRESSURE_vs_TEMPERATURE"]["rmse"], 1.0, places=12)

    def test_UnknownMetricRaises(self):
        with self.assertRaisesRegex(ValueError, "Unsupported metric"):
            self._CreateProcess(metrics='["mse", "psnr"]')

    def test_RelativeL2HandComputed(self):
        # predicted - reference = 1 per node; reference = (0,1,2,3)
        process = self._CreateProcess(metrics='["relative_l2"]')
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        expected = math.sqrt(4.0) / math.sqrt(0.0 + 1.0 + 4.0 + 9.0)
        self.assertAlmostEqual(
            process.history[0]["PRESSURE_vs_TEMPERATURE"]["relative_l2"], expected, places=12)

    def test_WeightedMetricsHandComputed(self):
        for node in self.model_part.Nodes:  # weights 1, 2, 3, 4
            node.SetValue(Kratos.NODAL_PAUX, float(node.Id))
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"     : "Main",
                "list_of_comparisons" : [
                    {
                        "predicted_variable" : "PRESSURE",
                        "predicted_location" : "node_historical",
                        "reference_variable" : "TEMPERATURE",
                        "reference_location" : "node_historical",
                        "weight_variable"    : "NODAL_PAUX",
                        "weight_location"    : "node_non_historical",
                        "metrics"            : ["weighted_mse", "weighted_rmse"]
                    }
                ],
                "output_file"         : "test_validation_metrics.json"
            }
        }""")
        process = validation_metrics_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        # per-node squared error = 1 everywhere -> weighted mean = 1 exactly
        values = process.history[0]["PRESSURE_vs_TEMPERATURE"]
        self.assertAlmostEqual(values["weighted_mse"], 1.0, places=12)
        self.assertAlmostEqual(values["weighted_rmse"], 1.0, places=12)

    def test_WeightedMetricWithoutWeightVariableRaises(self):
        with self.assertRaisesRegex(ValueError, "weight_variable"):
            self._CreateProcess(metrics='["weighted_mse"]')

    def test_ComputeMetricValuesWeighted(self):
        import torch
        predicted = torch.tensor([1.0, 2.0])
        reference = torch.tensor([0.0, 0.0])  # squared errors 1, 4
        weights = torch.tensor([3.0, 1.0])    # weighted mse = (3*1 + 1*4)/4
        values = validation_metrics_process.ComputeMetricValues(
            predicted, reference, ["weighted_mse"], weights)
        self.assertAlmostEqual(values["weighted_mse"], 7.0 / 4.0, places=12)
        with self.assertRaisesRegex(ValueError, "weights"):
            validation_metrics_process.ComputeMetricValues(
                predicted, reference, ["weighted_mse"])


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestEnsembleMetrics(KratosUnittest.TestCase):
    def test_CrpsOfTightAndLooseEnsembles(self):
        import torch
        torch.manual_seed(0)
        reference = torch.zeros(50)
        tight = 0.01 * torch.randn(8, 50)
        loose = 10.0 * torch.randn(8, 50) + 5.0

        tight_values = validation_metrics_process.ComputeEnsembleMetricValues(
            tight, reference, ["crps", "kcrps"])
        loose_values = validation_metrics_process.ComputeEnsembleMetricValues(
            loose, reference, ["crps", "kcrps"])
        for name in ("crps", "kcrps"):
            self.assertGreaterEqual(tight_values[name], 0.0)
            self.assertLess(tight_values[name], loose_values[name])

    def test_InvalidEnsembleShapesRejected(self):
        import torch
        with self.assertRaisesRegex(ValueError, "M >= 2"):
            validation_metrics_process.ComputeEnsembleMetricValues(
                torch.zeros(1, 5), torch.zeros(5), ["crps"])
        with self.assertRaisesRegex(ValueError, "Shape mismatch"):
            validation_metrics_process.ComputeEnsembleMetricValues(
                torch.zeros(3, 5), torch.zeros(4), ["crps"])
        with self.assertRaisesRegex(ValueError, "ensemble metric"):
            validation_metrics_process.ComputeEnsembleMetricValues(
                torch.zeros(3, 5), torch.zeros(5), ["spread"])


if __name__ == '__main__':
    KratosUnittest.main()
