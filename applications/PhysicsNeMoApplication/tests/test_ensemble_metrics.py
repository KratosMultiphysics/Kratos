"""Tests for the ensemble (CRPS) metrics settings block and for retaining the
ensemble an inference process would otherwise discard."""

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
from KratosMultiphysics.PhysicsNeMoApplication.processes import validation_metrics_process
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.metrics.general.crps import crps  # noqa: F401
    have_crps = have_torch
except ImportError:
    have_crps = False

_MISSING = "Missing required python modules: torch, physicsnemo."


@KratosUnittest.skipUnless(have_crps, _MISSING)
class TestEnsembleMetricValues(KratosUnittest.TestCase):

    def test_PerfectEnsembleScoresAboutZero(self):
        torch.manual_seed(0)
        reference = torch.randn(8, dtype=torch.float64)
        ensemble = reference.expand(16, 8).clone()
        values = validation_metrics_process.ComputeEnsembleMetricValues(
            ensemble, reference, ["crps", "kcrps"])
        # it can land marginally negative on float64 round-off
        self.assertLess(abs(values["crps"]), 1e-12)
        self.assertLess(abs(values["kcrps"]), 1e-12)

    def test_SpreadEnsembleScoresWorseThanATightOne(self):
        torch.manual_seed(0)
        reference = torch.zeros(64, dtype=torch.float64)
        tight = torch.randn(16, 64, dtype=torch.float64) * 0.1
        loose = torch.randn(16, 64, dtype=torch.float64) * 2.0
        tight_score = validation_metrics_process.ComputeEnsembleMetricValues(
            tight, reference, ["crps"])["crps"]
        loose_score = validation_metrics_process.ComputeEnsembleMetricValues(
            loose, reference, ["crps"])["crps"]
        self.assertLess(tight_score, loose_score)

    def test_Validation(self):
        reference = torch.zeros(4, dtype=torch.float64)
        with self.assertRaisesRegex(ValueError, "Unsupported ensemble metric"):
            validation_metrics_process.ComputeEnsembleMetricValues(
                torch.zeros(3, 4, dtype=torch.float64), reference, ["mse"])
        with self.assertRaisesRegex(ValueError, "M >= 2"):
            validation_metrics_process.ComputeEnsembleMetricValues(
                torch.zeros(1, 4, dtype=torch.float64), reference, ["crps"])


@KratosUnittest.skipUnless(have_crps, _MISSING)
class TestEnsembleComparisonsThroughProcess(KratosUnittest.TestCase):
    """The gap this closes: the metric shipped and was tested, but no
    settings key reached it, so it was unusable in a deployment."""

    def _Model(self, spread):
        model = Kratos.Model()
        model_part = model.CreateModelPart("Ens")
        model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        for i in range(32):
            node = model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.TEMPERATURE, 1.0)
            node.SetValue(Kratos.PRESSURE, 1.0 + spread)
            node.SetValue(Kratos.DISTANCE, 1.0 - spread)
            node.SetValue(Kratos.NODAL_ERROR, 1.0)
        model_part.ProcessInfo[Kratos.STEP] = 1
        model_part.ProcessInfo[Kratos.TIME] = 0.0
        return model, model_part

    def _Run(self, model, output_file):
        process = validation_metrics_process.Factory(Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"      : "Ens",
                "list_of_comparisons"  : [],
                "ensemble_comparisons" : [ {
                    "member_variables"   : ["PRESSURE", "DISTANCE", "NODAL_ERROR"],
                    "member_location"    : "node_non_historical",
                    "reference_variable" : "TEMPERATURE",
                    "metrics"            : ["crps", "kcrps"]
                } ],
                "output_file" : "%s"
            }
        }""" % output_file), model)
        try:
            process.ExecuteFinalizeSolutionStep()
            return process.history[-1]
        finally:
            KratosUtilities.DeleteFileIfExisting(output_file)

    def test_MetricReachesTheProcessAndRanksSpread(self):
        tight, _ = self._Model(0.05)
        loose, _ = self._Model(1.5)
        tight_record = self._Run(tight, "test_ens_tight.json")
        loose_record = self._Run(loose, "test_ens_loose.json")

        self.assertIn("ensemble_TEMPERATURE", tight_record)
        self.assertIn("crps", tight_record["ensemble_TEMPERATURE"])
        self.assertLess(tight_record["ensemble_TEMPERATURE"]["crps"],
                        loose_record["ensemble_TEMPERATURE"]["crps"])

    def test_Validation(self):
        model, _ = self._Model(0.1)
        with self.assertRaisesRegex(ValueError, "at least 2"):
            validation_metrics_process.Factory(Kratos.Parameters("""{
                "Parameters": {
                    "model_part_name"      : "Ens",
                    "list_of_comparisons"  : [],
                    "ensemble_comparisons" : [ { "member_variables" : ["PRESSURE"],
                                                 "reference_variable" : "TEMPERATURE" } ]
                }
            }"""), model)
        with self.assertRaisesRegex(ValueError, "Unsupported ensemble metric"):
            validation_metrics_process.Factory(Kratos.Parameters("""{
                "Parameters": {
                    "model_part_name"      : "Ens",
                    "list_of_comparisons"  : [],
                    "ensemble_comparisons" : [ { "member_variables" : ["PRESSURE", "DISTANCE"],
                                                 "reference_variable" : "TEMPERATURE",
                                                 "metrics" : ["rmse"] } ]
                }
            }"""), model)


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestRetainEnsemble(KratosUnittest.TestCase):
    """The (M, ...) stack exists in memory for one line before being reduced
    to mean/std; CRPS cannot be recovered from that reduction."""

    def setUp(self):
        self.checkpoints = [f"test_retain_member_{i}.pt" for i in range(3)]
        for index, path in enumerate(self.checkpoints):
            class Scaled(torch.nn.Module):
                def __init__(self, factor):
                    super().__init__()
                    self.factor = factor

                def forward(self, x):
                    return self.factor * x

            torch.jit.script(Scaled(1.0 + 0.5 * index)).save(path)

    def tearDown(self):
        for path in self.checkpoints:
            KratosUtilities.DeleteFileIfExisting(path)

    def _Process(self, model, retain):
        files = ", ".join(f'"{path}"' for path in self.checkpoints)
        return inference_process.Factory(Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Retain",
                "model_settings"  : { "checkpoint_files" : [%s],
                                      "checkpoint_type"  : "torchscript",
                                      "device"           : "cpu" },
                "input_fields"    : [ { "variable_name" : "PRESSURE",
                                        "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE",
                                        "data_location" : "node_historical" } ],
                "uncertainty"     : { "method" : "ensemble",
                                      "retain_ensemble" : %s,
                                      "uncertainty_fields" : [
                                          { "variable_name" : "NODAL_ERROR",
                                            "data_location" : "node_non_historical" } ] }
            }
        }""" % (files, "true" if retain else "false")), model)

    def _Model(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("Retain")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        for i in range(5):
            node = model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.PRESSURE, float(i) + 1.0)
        model_part.ProcessInfo[Kratos.STEP] = 1
        return model, model_part

    def test_MembersAreKeptWhenAsked(self):
        model, _ = self._Model()
        process = self._Process(model, retain=True)
        process.ExecuteFinalizeSolutionStep()
        self.assertIsNotNone(process.last_ensemble)
        self.assertEqual(tuple(process.last_ensemble.shape), (3, 5, 1))
        # the members really differ - the mean/std reduction loses this
        spread = process.last_ensemble.std(dim=0)
        self.assertGreater(float(spread.max()), 0.0)

    def test_NotKeptByDefault(self):
        model, _ = self._Model()
        process = self._Process(model, retain=False)
        process.ExecuteFinalizeSolutionStep()
        self.assertIsNone(process.last_ensemble)


if __name__ == '__main__':
    KratosUnittest.main()
