from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.active_learning import metrology

try:
    import physicsnemo.metrics.general.mse
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestValidationMetricsMetrology(KratosUnittest.TestCase):
    def setUp(self):
        self.output_file = Path("test_metrology_records.json")

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.output_file))

    def _CreateStrategy(self, pairs):
        settings = Kratos.Parameters("""{
            "metrics"     : ["mse", "rmse", "max_abs_error"],
            "output_file" : "%s"
        }""" % self.output_file)
        return metrology.CreateValidationMetricsMetrology(settings, lambda: pairs)

    def test_ComputeAppendsOneRecordPerCall(self):
        predicted = numpy.array([1.0, 2.0, 3.0])
        reference = numpy.array([1.0, 2.0, 5.0])
        strategy = self._CreateStrategy({"TEMPERATURE": (predicted, reference)})

        strategy.compute()
        strategy()  # __call__ dispatches to compute, as the driver does
        self.assertEqual(len(strategy.records), 2)
        self.assertEqual(len(strategy), 2)

        record = strategy.records[0]
        self.assertEqual(record["iteration"], 0)
        self.assertAlmostEqual(record["TEMPERATURE"]["mse"], 4.0 / 3.0, places=12)
        self.assertAlmostEqual(record["TEMPERATURE"]["rmse"], (4.0 / 3.0) ** 0.5, places=12)
        self.assertAlmostEqual(record["TEMPERATURE"]["max_abs_error"], 2.0, places=12)
        self.assertEqual(strategy.records[1]["iteration"], 1)

    def test_SerializeAndLoadRoundTrip(self):
        strategy = self._CreateStrategy({"T": (numpy.ones(4), numpy.zeros(4))})
        strategy.compute()
        strategy.serialize_records()
        self.assertTrue(self.output_file.is_file())

        restored = self._CreateStrategy({"T": (numpy.ones(4), numpy.zeros(4))})
        restored.load_records()
        self.assertEqual(restored.records, strategy.records)

        restored.reset()
        self.assertEqual(restored.records, [])

    def test_AttachProtocol(self):
        strategy = self._CreateStrategy({})
        self.assertFalse(strategy.is_attached)
        strategy.attach(object())
        self.assertTrue(strategy.is_attached)

    def test_UnsupportedMetricRaises(self):
        settings = Kratos.Parameters("""{ "metrics" : ["nonsense"] }""")
        with self.assertRaisesRegex(ValueError, "Unsupported metric"):
            metrology.CreateValidationMetricsMetrology(settings, lambda: {})

    def test_ShapeMismatchRaises(self):
        strategy = self._CreateStrategy({"T": (numpy.ones(3), numpy.ones(4))})
        with self.assertRaisesRegex(ValueError, "Shape mismatch"):
            strategy.compute()


if __name__ == '__main__':
    KratosUnittest.main()
