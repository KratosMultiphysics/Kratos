"""Tests for the GP-augmented uncertainty head and the calibration metrics
it feeds."""

import math

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
from KratosMultiphysics.PhysicsNeMoApplication.deployment import uncertainty_utils
from KratosMultiphysics.PhysicsNeMoApplication.processes import validation_metrics_process
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import physicsnemo.experimental.uq as _experimental_uq
    have_gp = have_torch and getattr(_experimental_uq, "_GPYTORCH_AVAILABLE", False)
except ImportError:
    have_gp = False

_MISSING = "Missing required python modules: torch, physicsnemo >= 2.2, gpytorch."


def _SyntheticFeatures(n=384, dim=4, seed=1):
    torch.manual_seed(seed)
    features = torch.randn(n, dim)
    targets = torch.stack([2.0 * features[:, 0], -1.0 * features[:, 1]], dim=1)
    return features, targets


@KratosUnittest.skipUnless(have_gp, _MISSING)
class TestGpHead(KratosUnittest.TestCase):

    def setUp(self):
        self.head_file = "test_gp_head.pt"
        self.features, self.targets = _SyntheticFeatures()
        self.settings = Kratos.Parameters("""{
            "num_tasks"      : 2,
            "n_inducing"     : 48,
            "mlp_hidden"     : [16],
            "epochs"         : 150,
            "learning_rate"  : 0.05,
            "kl_ramp_epochs" : 40,
            "seed"           : 0
        }""")

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(self.head_file)

    def test_FitReducesTheElboAndTracksTheFunction(self):
        head, history = uncertainty_utils.FitGpHead(
            self.features, self.targets, self.settings)
        self.assertEqual(len(history), 150)
        self.assertLess(history[-1], history[0])

        mean, std = uncertainty_utils.PredictWithGpHead(head, self.features)
        self.assertEqual(tuple(mean.shape), (384, 2))
        self.assertEqual(mean.dtype, torch.float64)
        self.assertLess(float(((mean - self.targets.double()) ** 2).mean()), 0.05)
        self.assertTrue(bool((std > 0).all()))

    def test_EpistemicVarianceGrowsAwayFromTheTrainingData(self):
        """The property that makes the head worth having: an error metric
        cannot tell you where the model is guessing, but this can."""
        head, _ = uncertainty_utils.FitGpHead(
            self.features, self.targets, self.settings)

        _, inside = uncertainty_utils.PredictWithGpHead(
            head, self.features, epistemic_only=True)
        torch.manual_seed(7)
        far_away = torch.randn(96, 4) * 6.0
        _, outside = uncertainty_utils.PredictWithGpHead(
            head, far_away, epistemic_only=True)
        self.assertGreater(float(outside.mean()), float(inside.mean()))

    def test_SaveAndLoadRoundTrip(self):
        head, _ = uncertainty_utils.FitGpHead(
            self.features, self.targets, self.settings)
        expected_mean, expected_std = uncertainty_utils.PredictWithGpHead(
            head, self.features)

        uncertainty_utils.SaveGpHead(head, self.head_file, config={
            "input_dim": 4, "n_train": 384, "num_tasks": 2,
            "n_inducing": 48, "mlp_hidden": [16]})
        restored = uncertainty_utils.LoadGpHead(self.head_file)
        mean, std = uncertainty_utils.PredictWithGpHead(restored, self.features)

        numpy.testing.assert_allclose(mean.numpy(), expected_mean.numpy(), atol=1e-10)
        numpy.testing.assert_allclose(std.numpy(), expected_std.numpy(), atol=1e-10)

    def test_Validation(self):
        with self.assertRaisesRegex(ValueError, "point count"):
            uncertainty_utils.FitGpHead(self.features, self.targets[:10], self.settings)
        with self.assertRaisesRegex(ValueError, r"\(N, D\)"):
            uncertainty_utils.FitGpHead(self.features[0], self.targets, self.settings)
        KratosUtilities.DeleteFileIfExisting("no_config.pt")
        try:
            head, _ = uncertainty_utils.FitGpHead(
                self.features, self.targets, self.settings)
            uncertainty_utils.SaveGpHead(head, "no_config.pt")
            with self.assertRaisesRegex(ValueError, "input_dim"):
                uncertainty_utils.LoadGpHead("no_config.pt")
        finally:
            KratosUtilities.DeleteFileIfExisting("no_config.pt")


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestCalibrationMetrics(KratosUnittest.TestCase):
    """Analytic pins: a perfectly calibrated Gaussian has known values."""

    def test_CoverageAndNllMatchTheClosedForm(self):
        torch.manual_seed(0)
        n = 200000
        sigma = 0.5
        mean = torch.zeros(n, dtype=torch.float64)
        std = torch.full((n,), sigma, dtype=torch.float64)
        reference = torch.randn(n, dtype=torch.float64) * sigma

        values = validation_metrics_process.ComputeCalibrationMetricValues(
            mean, std, reference, ["coverage", "nll", "sharpness", "calibration_error"])

        self.assertAlmostEqual(values["coverage"], 0.95, delta=0.01)
        # closed form for a correctly specified Gaussian: 0.5*(1+ln(2*pi*s^2))
        expected_nll = 0.5 * (1.0 + math.log(2.0 * math.pi * sigma ** 2))
        self.assertAlmostEqual(values["nll"], expected_nll, delta=0.02)
        self.assertAlmostEqual(values["sharpness"], sigma, places=12)
        self.assertLess(values["calibration_error"], 0.01)

    def test_OverconfidenceIsPunished(self):
        torch.manual_seed(0)
        n = 50000
        reference = torch.randn(n, dtype=torch.float64)
        mean = torch.zeros(n, dtype=torch.float64)
        honest = validation_metrics_process.ComputeCalibrationMetricValues(
            mean, torch.ones(n, dtype=torch.float64), reference, ["coverage", "nll"])
        overconfident = validation_metrics_process.ComputeCalibrationMetricValues(
            mean, torch.full((n,), 0.2, dtype=torch.float64), reference, ["coverage", "nll"])
        self.assertLess(overconfident["coverage"], honest["coverage"])
        self.assertGreater(overconfident["nll"], honest["nll"])

    def test_Validation(self):
        values = torch.zeros(4, dtype=torch.float64)
        with self.assertRaisesRegex(ValueError, "Unsupported calibration metric"):
            validation_metrics_process.ComputeCalibrationMetricValues(
                values, values + 1.0, values, ["bogus"])
        with self.assertRaisesRegex(ValueError, "Shape mismatch"):
            validation_metrics_process.ComputeCalibrationMetricValues(
                values, values + 1.0, torch.zeros(3, dtype=torch.float64), ["coverage"])
        with self.assertRaisesRegex(ValueError, "negative"):
            validation_metrics_process.ComputeCalibrationMetricValues(
                values, values - 1.0, values, ["coverage"])


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestCalibrationThroughValidationProcess(KratosUnittest.TestCase):

    def test_UncertaintyComparisonsReachTheMetrics(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("Calib")
        model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        for i in range(64):
            node = model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.TEMPERATURE, 1.0)
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0)
            node.SetValue(Kratos.NODAL_ERROR, 0.5)

        process = validation_metrics_process.Factory(Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Calib",
                "list_of_comparisons" : [],
                "uncertainty_comparisons" : [ {
                    "mean_variable"      : "TEMPERATURE",
                    "std_variable"       : "NODAL_ERROR",
                    "std_location"       : "node_non_historical",
                    "reference_variable" : "PRESSURE",
                    "metrics"            : ["coverage", "sharpness"]
                } ],
                "output_file" : "test_calibration_metrics.json"
            }
        }"""), model)
        model_part.ProcessInfo[Kratos.STEP] = 1
        model_part.ProcessInfo[Kratos.TIME] = 0.0
        try:
            process.ExecuteFinalizeSolutionStep()
            record = process.history[-1]
            key = "calibration_TEMPERATURE_vs_PRESSURE"
            self.assertIn(key, record)
            # mean == reference exactly, so every point is covered
            self.assertAlmostEqual(record[key]["coverage"], 1.0, places=12)
            self.assertAlmostEqual(record[key]["sharpness"], 0.5, places=12)
        finally:
            KratosUtilities.DeleteFileIfExisting("test_calibration_metrics.json")

    def test_UnsupportedCalibrationMetricRejected(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("CalibBad")
        model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        with self.assertRaisesRegex(ValueError, "Unsupported calibration metric"):
            validation_metrics_process.Factory(Kratos.Parameters("""{
                "Parameters": {
                    "model_part_name" : "CalibBad",
                    "list_of_comparisons" : [],
                    "uncertainty_comparisons" : [ {
                        "mean_variable"      : "TEMPERATURE",
                        "std_variable"       : "TEMPERATURE",
                        "reference_variable" : "TEMPERATURE",
                        "metrics"            : ["rmse"]
                    } ]
                }
            }"""), model)


@KratosUnittest.skipUnless(have_gp, _MISSING)
class TestGpUncertaintyThroughProcess(KratosUnittest.TestCase):

    def setUp(self):
        self.checkpoint = "test_gp_surrogate.pt"
        self.head_file = "test_gp_surrogate.pt.gp_head.pt"

    def tearDown(self):
        for path in (self.checkpoint, self.head_file, self.checkpoint + ".card.json"):
            KratosUtilities.DeleteFileIfExisting(path)

    def test_GpMethodWritesPerNodeUncertainty(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("GpDeploy")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        for i in range(48):
            node = model_part.CreateNewNode(i + 1, float(i) / 48.0, 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.PRESSURE, float(i) / 48.0)

        torch.manual_seed(0)
        # the gathered fields are float64, so the deployed model must be too
        surrogate = torch.nn.Sequential(torch.nn.Linear(1, 8), torch.nn.Tanh(),
                                        torch.nn.Linear(8, 1)).double()
        training_utils.SaveTrainedModel(surrogate, self.checkpoint)

        features = torch.rand(256, 1)
        head, _ = uncertainty_utils.FitGpHead(
            features, features * 2.0, Kratos.Parameters("""{
                "num_tasks": 1, "n_inducing": 32, "mlp_hidden": [16],
                "epochs": 60, "seed": 0
            }"""))
        uncertainty_utils.SaveGpHead(head, self.head_file, config={
            "input_dim": 1, "n_train": 256, "num_tasks": 1,
            "n_inducing": 32, "mlp_hidden": [16]})

        process = inference_process.Factory(Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "GpDeploy",
                "model_settings"  : { "checkpoint_file" : "%s",
                                      "checkpoint_type" : "torchscript",
                                      "device"          : "cpu" },
                "input_fields"    : [ { "variable_name" : "PRESSURE",
                                        "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE",
                                        "data_location" : "node_historical" } ],
                "uncertainty"     : {
                    "method"             : "gp",
                    "gp_head_file"       : "%s",
                    "uncertainty_fields" : [ { "variable_name" : "NODAL_ERROR",
                                               "data_location" : "node_non_historical" } ]
                }
            }
        }""" % (self.checkpoint, self.head_file)), model)
        model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        spread = numpy.array([n.GetValue(Kratos.NODAL_ERROR) for n in model_part.Nodes])
        self.assertEqual(spread.shape, (48,))
        self.assertTrue(numpy.isfinite(spread).all())
        self.assertTrue((spread > 0.0).all())

    def test_GpMethodNeedsAHeadFile(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("GpNoHead")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        with self.assertRaisesRegex(ValueError, "gp_head_file"):
            inference_process.Factory(Kratos.Parameters("""{
                "Parameters": {
                    "model_part_name" : "GpNoHead",
                    "model_settings"  : { "checkpoint_file" : "unused.pt" },
                    "input_fields"    : [ { "variable_name" : "PRESSURE",
                                            "data_location" : "node_historical" } ],
                    "output_fields"   : [ { "variable_name" : "PRESSURE",
                                            "data_location" : "node_historical" } ],
                    "uncertainty"     : { "method" : "gp" }
                }
            }"""), model)


if __name__ == '__main__':
    KratosUnittest.main()
