"""Tests for the validation metrics that look past pointwise error: the
distribution metrics (histogram L1, entropy difference, relative MSE) and
the spectral comparisons a superresolution model should be judged on."""

import math

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.processes import validation_metrics_process

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.metrics.general.power_spectrum import power_spectrum  # noqa: F401
    from physicsnemo.metrics.general.histogram import histogram  # noqa: F401
    have_metrics = True
except ImportError:
    have_metrics = False


def _SmoothPlane(size=32):
    y, x = torch.meshgrid(torch.linspace(0.0, 1.0, size), torch.linspace(0.0, 1.0, size),
                          indexing="ij")
    return torch.sin(2.0 * math.pi * x) * torch.cos(2.0 * math.pi * y)


@KratosUnittest.skipUnless(have_torch and have_metrics,
                           "Missing required python modules: torch, physicsnemo.")
class TestDistributionMetrics(KratosUnittest.TestCase):
    def test_IdenticalFieldsScoreZero(self):
        torch.manual_seed(0)
        field = torch.randn(500, 1, dtype=torch.float64)
        values = validation_metrics_process.ComputeMetricValues(
            field, field.clone(), ["relative_mse", "histogram_l1", "entropy_difference"])
        for name, value in values.items():
            self.assertAlmostEqual(value, 0.0, places=12, msg=name)

    def test_RelativeMseIsTheSquareOfRelativeL2(self):
        torch.manual_seed(1)
        reference = torch.randn(200, 3, dtype=torch.float64)
        predicted = reference + 0.1 * torch.randn(200, 3, dtype=torch.float64)
        values = validation_metrics_process.ComputeMetricValues(
            predicted, reference, ["relative_mse", "relative_l2"])
        self.assertAlmostEqual(values["relative_mse"], values["relative_l2"] ** 2, places=12)

    def test_DisjointDistributionsAreAtTheL1Maximum(self):
        predicted = torch.zeros(400, 1, dtype=torch.float64)
        reference = torch.ones(400, 1, dtype=torch.float64)
        values = validation_metrics_process.ComputeMetricValues(
            predicted, reference, ["histogram_l1"], bins=16)
        self.assertAlmostEqual(values["histogram_l1"], 2.0, places=10)

    def test_ACollapsedFieldIsCaughtWhereTheMeanIsNot(self):
        """What these metrics are for: a prediction that is right on
        average but has lost the field's extremes."""
        torch.manual_seed(2)
        reference = torch.randn(2000, 1, dtype=torch.float64)
        collapsed = 0.2 * reference  # same mean, a fifth of the spread
        values = validation_metrics_process.ComputeMetricValues(
            collapsed, reference, ["histogram_l1", "entropy_difference"], bins=32)
        self.assertAlmostEqual(float(collapsed.mean()), float(0.2 * reference.mean()), places=12)
        self.assertGreater(values["histogram_l1"], 0.5)
        self.assertGreater(values["entropy_difference"], 0.5)

    def test_UpstreamsNormalizedEntropyLeavesItsDocumentedRange(self):
        """Why the metric uses UNNORMALIZED entropy: the normalized form is
        documented to map onto [0, 1] and returns a negative value for a
        delta distribution."""
        from physicsnemo.metrics.general.entropy import entropy_from_counts

        delta = torch.zeros(20, dtype=torch.float64)
        delta[3] = 1.0
        edges = torch.linspace(0.0, 1.0, 21, dtype=torch.float64)
        self.assertLess(float(entropy_from_counts(delta, edges, normalized=True)), 0.0)

    def test_UnknownMetricRaises(self):
        field = torch.zeros(4, 1)
        with self.assertRaisesRegex(ValueError, "Unsupported metric"):
            validation_metrics_process.ComputeMetricValues(field, field, ["kl_divergence"])


@KratosUnittest.skipUnless(have_torch and have_metrics,
                           "Missing required python modules: torch, physicsnemo.")
class TestSpectralMetrics(KratosUnittest.TestCase):
    def test_IdenticalGridsHaveNoSpectralError(self):
        grid = _SmoothPlane()[None].to(torch.float64)
        values, spectra = validation_metrics_process.ComputeSpectralMetricValues(
            grid, grid.clone(), list(validation_metrics_process.SPECTRAL_METRICS))
        self.assertAlmostEqual(values["power_spectrum_relative_l2"], 0.0, places=12)
        self.assertAlmostEqual(values["high_wavenumber_energy_ratio"], 1.0, places=12)
        self.assertEqual(len(spectra["wavenumbers"]), len(spectra["predicted"]))
        self.assertEqual(len(spectra["predicted"]), len(spectra["reference"]))

    def test_ASmoothPredictionOfARoughFieldIsTooSmooth(self):
        """The superresolution failure pointwise error hides: the smooth
        answer is close everywhere and carries almost none of the reference's
        fine-scale energy."""
        torch.manual_seed(3)
        smooth = _SmoothPlane().to(torch.float64)
        reference = (smooth + 0.15 * torch.randn_like(smooth))[None]
        predicted = smooth[None]

        pointwise = validation_metrics_process.ComputeMetricValues(
            predicted, reference, ["rmse"])
        values, _ = validation_metrics_process.ComputeSpectralMetricValues(
            predicted, reference, ["high_wavenumber_energy_ratio"])
        self.assertLess(pointwise["rmse"], 0.2)                    # looks fine pointwise
        self.assertLess(values["high_wavenumber_energy_ratio"], 0.1)  # but is far too smooth

    def test_LeadingAxesAreAveraged(self):
        grid = _SmoothPlane().to(torch.float64)
        stacked = torch.stack([grid, grid, grid])  # (C, H, W)
        single, _ = validation_metrics_process.ComputeSpectralMetricValues(
            grid[None], grid[None], ["power_spectrum_relative_l2"])
        channels, _ = validation_metrics_process.ComputeSpectralMetricValues(
            stacked, stacked, ["power_spectrum_relative_l2"])
        self.assertAlmostEqual(single["power_spectrum_relative_l2"],
                               channels["power_spectrum_relative_l2"], places=12)

    def test_Validation(self):
        grid = _SmoothPlane()[None]
        with self.assertRaisesRegex(ValueError, "Shape mismatch"):
            validation_metrics_process.ComputeSpectralMetricValues(
                grid, grid[..., :16], ["power_spectrum_relative_l2"])
        with self.assertRaisesRegex(ValueError, "4x4"):
            validation_metrics_process.ComputeSpectralMetricValues(
                torch.zeros(1, 2, 2), torch.zeros(1, 2, 2), ["power_spectrum_relative_l2"])
        with self.assertRaisesRegex(ValueError, "spectral metric"):
            validation_metrics_process.ComputeSpectralMetricValues(
                grid, grid, ["phase_error"])


@KratosUnittest.skipUnless(have_torch and have_metrics,
                           "Missing required python modules: torch, physicsnemo.")
class TestSpectralComparisonsThroughTheProcess(KratosUnittest.TestCase):
    def setUp(self):
        from test_grid_bridge import CreateStructuredTetModelPart

        self.report = "test_spectral_report.json"
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=8,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        rng = numpy.random.default_rng(0)
        for node in self.model_part.Nodes:
            smooth = math.sin(2.0 * math.pi * node.X) * math.cos(2.0 * math.pi * node.Y)
            node.SetSolutionStepValue(Kratos.TEMPERATURE, smooth)                      # predicted
            node.SetSolutionStepValue(Kratos.PRESSURE, smooth + 0.3 * rng.standard_normal())  # reference

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting(self.report)

    def _Process(self, spectral_block, comparisons_block=None):
        comparisons = comparisons_block or """[
            {
                "predicted_variable" : "TEMPERATURE",
                "reference_variable" : "PRESSURE",
                "metrics"            : ["rmse"]
            }
        ]"""
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"      : "Main",
                "list_of_comparisons"  : %s,
                "spectral_comparisons" : %s,
                "output_file"          : "%s"
            }
        }""" % (comparisons, spectral_block, self.report))
        return validation_metrics_process.Factory(settings, self.model)

    def test_TheSpectrumIsRecordedAndSaysTooSmooth(self):
        process = self._Process("""[
            {
                "predicted_variable" : "TEMPERATURE",
                "reference_variable" : "PRESSURE",
                "grid_shape"         : [16, 16, 2],
                "squeeze_axis"       : 2
            }
        ]""")
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        record = process.history[-1]
        self.assertIn("spectrum_TEMPERATURE_vs_PRESSURE", record)
        entry = record["spectrum_TEMPERATURE_vs_PRESSURE"]
        self.assertLess(entry["high_wavenumber_energy_ratio"], 1.0)
        self.assertGreater(len(entry["spectra"]["wavenumbers"]), 0)

    def test_AnUnconfiguredProcessHasNoPhantomSpectralEntry(self):
        """The list-defaults trap: a one-entry schema as the default would
        be inherited by every existing configuration as a real entry."""
        process = self._Process("[]")
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        self.assertFalse(any(key.startswith("spectrum_") for key in process.history[-1]))

    def test_BinsReachTheDistributionMetrics(self):
        process = self._Process("[]", """[
            {
                "predicted_variable" : "TEMPERATURE",
                "reference_variable" : "PRESSURE",
                "bins"               : 8,
                "metrics"            : ["histogram_l1", "entropy_difference", "relative_mse"]
            }
        ]""")
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        values = process.history[-1]["TEMPERATURE_vs_PRESSURE"]
        self.assertGreater(values["histogram_l1"], 0.0)
        self.assertGreater(values["relative_mse"], 0.0)

    def test_InvalidSpectralSettingsRaise(self):
        with self.assertRaisesRegex(ValueError, "spectral metric"):
            self._Process("""[ { "predicted_variable" : "TEMPERATURE",
                                 "reference_variable" : "PRESSURE",
                                 "metrics" : ["phase_error"] } ]""")
        with self.assertRaisesRegex(ValueError, "grid_shape"):
            self._Process("""[ { "predicted_variable" : "TEMPERATURE",
                                 "reference_variable" : "PRESSURE",
                                 "grid_shape" : [16, 16] } ]""")


if __name__ == '__main__':
    KratosUnittest.main()
