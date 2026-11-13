from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.deployment import ood_guard_utils
from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import physicsnemo.experimental.guardrails.embedded  # noqa: F401
    have_guardrails = True
except ImportError:
    have_guardrails = False


def _MakeCalibratedGuard(n_samples=32, n_points=20, width=2, seed=0):
    torch.manual_seed(seed)
    samples = [torch.randn(n_points, width) for _ in range(n_samples)]
    guard = ood_guard_utils.CreateOODGuard(buffer_size=n_samples, feature_width=width)
    ood_guard_utils.CalibrateGuardFromTensors(guard, samples)
    return guard


@KratosUnittest.skipUnless(have_torch and have_guardrails,
                           "Missing required python modules: torch, physicsnemo (experimental).")
class TestOODGuardUtils(KratosUnittest.TestCase):
    def test_InDistributionPasses(self):
        guard = _MakeCalibratedGuard()
        torch.manual_seed(1)
        messages = ood_guard_utils.CheckFeatures(guard, torch.randn(20, 2) * 0.5)
        self.assertEqual(messages, [])

    def test_OutOfDistributionFlagged(self):
        guard = _MakeCalibratedGuard()
        torch.manual_seed(1)
        messages = ood_guard_utils.CheckFeatures(guard, torch.randn(20, 2) * 100.0)
        self.assertGreater(len(messages), 0)

    def test_SaveLoadRoundTrip(self):
        guard = _MakeCalibratedGuard()
        path = Path("test_ood_guard_roundtrip.pt")
        try:
            ood_guard_utils.SaveGuard(guard, path)
            loaded = ood_guard_utils.LoadGuard(path)
            torch.manual_seed(1)
            in_dist = torch.randn(20, 2) * 0.5
            out_dist = torch.randn(20, 2) * 100.0
            self.assertEqual(ood_guard_utils.CheckFeatures(loaded, in_dist), [])
            self.assertGreater(len(ood_guard_utils.CheckFeatures(loaded, out_dist)), 0)
        finally:
            KratosUtilities.DeleteFileIfExisting(str(path))

    def test_GuardCheckPolicies(self):
        guard = _MakeCalibratedGuard()
        path = Path("test_ood_guard_policies.pt")
        try:
            ood_guard_utils.SaveGuard(guard, path)
            torch.manual_seed(1)
            out_dist = torch.randn(20, 2) * 100.0

            advisory = ood_guard_utils.GuardCheck(Kratos.Parameters(
                '{"guard_file": "%s"}' % path))
            self.assertTrue(advisory.Check(out_dist, "Test"))  # warns, no raise
            self.assertTrue(advisory.last_flagged)

            strict = ood_guard_utils.GuardCheck(Kratos.Parameters(
                '{"guard_file": "%s", "policy": "strict"}' % path))
            with self.assertRaisesRegex(RuntimeError, "out-of-distribution"):
                strict.Check(out_dist, "Test")

            ignore = ood_guard_utils.GuardCheck(Kratos.Parameters(
                '{"guard_file": "%s", "policy": "ignore"}' % path))
            self.assertFalse(ignore.Check(out_dist, "Test"))
            self.assertFalse(ignore.enabled)

            disabled = ood_guard_utils.GuardCheck(Kratos.Parameters("{}"))
            self.assertFalse(disabled.Check(out_dist, "Test"))
        finally:
            KratosUtilities.DeleteFileIfExisting(str(path))

    def test_UnknownPolicyRaises(self):
        with self.assertRaisesRegex(ValueError, "policy"):
            ood_guard_utils.GuardCheck(Kratos.Parameters('{"policy": "loose"}'))


@KratosUnittest.skipUnless(have_torch and have_guardrails,
                           "Missing required python modules: torch, physicsnemo (experimental).")
class TestOODGuardEndToEnd(KratosUnittest.TestCase):
    """TrainModel calibrates + saves the guard; InferenceProcess checks it."""

    def setUp(self):
        self.checkpoint = Path("test_ood_e2e_model.pt")
        self.guard_file = Path("test_ood_e2e_model.pt.ood_guard.pt")

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))
        KratosUtilities.DeleteFileIfExisting(str(self.guard_file))

    def test_GridCalibrationMatchesTheGridCheckLayout(self):
        """Grid batches are (B, C, *spatial), channels FIRST.

        Calibration used to pool over everything but the LAST axis - for a
        grid that is a spatial axis, not the channels - so the guard scored
        a different quantity than the grid processes later check and flagged
        every input, in-distribution ones included.
        """
        torch.manual_seed(0)
        # channels with DISJOINT ranges (10..11 and 0..1): under the wrong
        # layout the rows pair same-channel values, so an in-distribution
        # (c0, c1) ~ (10.5, 0.5) point is far from both same-channel
        # clusters and gets flagged - identical uniform channels would make
        # both layouts statistically indistinguishable and prove nothing
        inputs = torch.rand(24, 2, 8, 8, 2)
        inputs[:, 0] += 10.0
        dataset = torch.utils.data.TensorDataset(inputs, torch.rand(24, 1, 8, 8, 2))
        # batch_size 8 gives three collect calls: with a single call the
        # guard's kNN check is inert and this test cannot discriminate the
        # layouts (verified by probing both) - the mutation check caught a
        # first version of this test being exactly that vacuous
        training_utils.TrainModel(torch.nn.Conv3d(2, 1, 1), dataset, Kratos.Parameters("""{
            "epochs"     : 1,
            "seed"       : 0,
            "batch_size" : 8,
            "device"     : "cpu",
            "ood_guard"  : { "guard_file" : "%s", "sensitivity" : 6.0 }
        }""" % self.guard_file))
        guard = ood_guard_utils.LoadGuard(str(self.guard_file))

        # what the grid deployment processes feed Check: (prod(spatial), C)
        def AsCheckFeatures(grid):
            return grid.reshape(grid.shape[0], -1).T.to(torch.float32)

        probe = torch.rand(2, 8, 8, 2)
        probe[0] += 10.0
        in_distribution = AsCheckFeatures(probe)
        far_out = AsCheckFeatures(probe + 100.0)
        self.assertFalse(ood_guard_utils.CheckFeatures(guard, in_distribution),
                         "an in-distribution grid was flagged - the calibration "
                         "layout does not match the check layout")
        self.assertTrue(ood_guard_utils.CheckFeatures(guard, far_out))

    def test_GuardCalibrationOnCudaTraining(self):
        """The guard is a CPU deployment artifact; training on CUDA must
        still calibrate it. collect() used to receive device tensors and
        died with "Expected all tensors to be on the same device"."""
        if not torch.cuda.is_available():
            self.skipTest("Requires a CUDA device.")

        torch.manual_seed(0)
        inputs = torch.rand(16, 1)
        dataset = torch.utils.data.TensorDataset(inputs, 2.0 * inputs)
        training_utils.TrainModel(torch.nn.Linear(1, 1), dataset, Kratos.Parameters("""{
            "epochs"    : 2,
            "seed"      : 0,
            "device"    : "cuda",
            "ood_guard" : { "guard_file" : "%s" }
        }""" % self.guard_file))
        self.assertTrue(self.guard_file.is_file())

    def test_TrainCalibratesGuardAndProcessChecksIt(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
        torch.manual_seed(0)
        inputs = torch.rand(64, 1, dtype=torch.float64)  # training inputs in [0, 1]
        dataset = torch.utils.data.TensorDataset(inputs, 2.0 * inputs)
        model = torch.nn.Linear(1, 1).double()
        training_utils.TrainModel(model, dataset, Kratos.Parameters("""{
            "epochs"    : 2,
            "seed"      : 0,
            "device"    : "cpu",
            "ood_guard" : { "guard_file" : "%s" }
        }""" % self.guard_file))
        self.assertTrue(self.guard_file.is_file())
        torch.jit.script(model).save(str(self.checkpoint))

        kratos_model = Kratos.Model()
        model_part = CreateStructuredTetModelPart(
            kratos_model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        settings_template = """{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : { "checkpoint_file" : "%s", "device" : "cpu" },
                "ood_guard"       : { "guard_file" : "%s", "policy" : "%s" },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }"""

        # in-distribution inputs: advisory guard stays silent
        for node in model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 0.5)
        process = inference_process.Factory(Kratos.Parameters(
            settings_template % (self.checkpoint, self.guard_file, "advisory")), kratos_model)
        model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        self.assertFalse(process._ood_guard.last_flagged)

        # far-out-of-distribution inputs: strict guard refuses to run
        for node in model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0e4)
        process = inference_process.Factory(Kratos.Parameters(
            settings_template % (self.checkpoint, self.guard_file, "strict")), kratos_model)
        with self.assertRaisesRegex(RuntimeError, "out-of-distribution"):
            process.ExecuteFinalizeSolutionStep()


if __name__ == '__main__':
    KratosUnittest.main()
