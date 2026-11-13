from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.deployment import onnx_utils
from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.deploy.onnx import export_to_onnx_stream  # noqa: F401
    have_onnx_export = True
except ImportError:
    have_onnx_export = False

try:
    import onnxruntime  # noqa: F401
    have_onnxruntime = True
except ImportError:
    have_onnxruntime = False

# A CUDA device is not enough: the installed onnxruntime must be a CUDA build.
# The reference environment deliberately keeps the CPU build (onnxruntime and
# onnxruntime-gpu install the same package directory and would clobber each
# other), so these self-skip there.
have_ort_cuda = (have_onnxruntime
                 and "CUDAExecutionProvider" in onnxruntime.get_available_providers())


class TestExportOnnxModelValidation(KratosUnittest.TestCase):
    def test_ExtensionEnforced(self):
        # the extension check runs before any optional import
        with self.assertRaisesRegex(ValueError, r"\.onnx"):
            training_utils.ExportOnnxModel(None, None, "model.pt")


class TestOnnxDeviceParsing(KratosUnittest.TestCase):
    """Device-string handling. Needs no GPU and no onnxruntime, so it runs
    everywhere - which matters, because the bug it guards (a dropped device
    index) is invisible on a single-GPU machine."""

    def test_ParsesCpuAndCuda(self):
        self.assertEqual(onnx_utils.ParseDevice("cpu"), (False, 0))
        self.assertEqual(onnx_utils.ParseDevice("cuda"), (True, 0))
        self.assertEqual(onnx_utils.ParseDevice("cuda:0"), (True, 0))
        self.assertEqual(onnx_utils.ParseDevice("cuda:3"), (True, 3))
        self.assertEqual(onnx_utils.ParseDevice("CUDA:2"), (True, 2))

    def test_MalformedDeviceRejected(self):
        # previously "cuda:x" was accepted by a substring test and the index
        # silently discarded, so the session ran on device 0
        for bad in ("cuda:", "cuda:x", "cuda-1", "cuda:1:2"):
            with self.subTest(device=bad):
                with self.assertRaisesRegex(ValueError, "cuda:N"):
                    onnx_utils.ParseDevice(bad)


@KratosUnittest.skipUnless(have_torch and have_onnx_export and have_onnxruntime,
                           "Missing required python modules: torch, physicsnemo, onnxruntime.")
class TestOnnxExportAndSession(KratosUnittest.TestCase):
    def setUp(self):
        self.onnx_file = Path("test_onnx_export_model.onnx")

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.onnx_file))
        KratosUtilities.DeleteFileIfExisting(str(self.onnx_file) + ".card.json")

    def test_ExportedModelMatchesTorch(self):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import onnx_utils
        torch.manual_seed(0)
        model = torch.nn.Sequential(
            torch.nn.Linear(4, 8), torch.nn.Tanh(), torch.nn.Linear(8, 2))
        sample = torch.randn(5, 4)
        training_utils.ExportOnnxModel(model, sample, self.onnx_file,
                                       card={"notes": "onnx export test"})
        self.assertTrue(self.onnx_file.is_file())
        self.assertTrue(Path(str(self.onnx_file) + ".card.json").is_file())

        with torch.no_grad():
            expected = model(sample).numpy()
        session = onnx_utils.CreateOrtSession(self.onnx_file, "cpu")
        (actual,) = session.run(None, {session.get_inputs()[0].name: sample.numpy()})
        self.assertLess(numpy.abs(expected - actual).max(), 1e-5)


@KratosUnittest.skipUnless(have_torch and have_onnx_export and have_onnxruntime,
                           "Missing required python modules: torch, physicsnemo, onnxruntime.")
class TestOnnxInferenceProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.onnx_file = Path("test_onnx_inference_model.onnx")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 10.0 * node.X + node.Y)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.onnx_file))
        KratosUtilities.DeleteFileIfExisting(str(self.onnx_file) + ".card.json")

    def _ExportAffineModel(self, card=None):
        # TEMPERATURE = 2 * PRESSURE + 1, exactly checkable through ORT
        model = torch.nn.Linear(1, 1)
        with torch.no_grad():
            model.weight.fill_(2.0)
            model.bias.fill_(1.0)
        n_nodes = self.model_part.NumberOfNodes()
        training_utils.ExportOnnxModel(
            model, torch.zeros(n_nodes, 1), self.onnx_file, card=card)

    def _CreateProcess(self, model_card_policy="advisory"):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import onnx_inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {
                    "onnx_file"         : "%s",
                    "device"            : "cpu",
                    "model_card_policy" : "%s"
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""" % (self.onnx_file, model_card_policy))
        return onnx_inference_process.Factory(settings, self.model)

    def test_ProcessRoundTrip(self):
        self._ExportAffineModel()
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE),
                2.0 * (10.0 * node.X + node.Y) + 1.0, places=5)

    def test_SessionIsCached(self):
        self._ExportAffineModel()
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        session = process._session
        self.assertIsNotNone(session)
        process.ExecuteFinalizeSolutionStep()
        self.assertIs(process._session, session)

    def test_MismatchedCardIsAdvisory(self):
        self._ExportAffineModel(card={
            "input_fields": [{"variable_name": "VELOCITY", "data_location": "node_historical"}]})
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()  # advisory: must not raise
        self.assertNotEqual(
            self.model_part.GetNode(1).GetSolutionStepValue(Kratos.TEMPERATURE), 0.0)

    def test_MismatchedCardStrictRaises(self):
        self._ExportAffineModel(card={
            "input_fields": [{"variable_name": "VELOCITY", "data_location": "node_historical"}]})
        process = self._CreateProcess(model_card_policy="strict")
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        with self.assertRaisesRegex(RuntimeError, "strict"):
            process.ExecuteFinalizeSolutionStep()

    def test_OutputNormalizationFromTheCardIsApplied(self):
        # this process keys its card off "onnx_file", not "checkpoint_file";
        # the override that makes the lookup work had no test
        mean, std = 30.0, 5.0
        self._ExportAffineModel(card={
            "output_normalization": {"type": "mean_std", "mean": [mean], "std": [std]}})
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE),
                std * (2.0 * (10.0 * node.X + node.Y) + 1.0) + mean, places=4)

    def test_InputNormalizationFromTheCardIsApplied(self):
        in_mean, in_std = 5.0, 2.0
        self._ExportAffineModel(card={
            "input_normalization": {"type": "mean_std", "mean": [in_mean], "std": [in_std]}})
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE),
                2.0 * ((10.0 * node.X + node.Y) - in_mean) / in_std + 1.0, places=4)

    def test_MultiOutputModelRaises(self):
        class TwoOutputs(torch.nn.Module):
            def forward(self, x):
                return x, 2.0 * x

        n_nodes = self.model_part.NumberOfNodes()
        training_utils.ExportOnnxModel(
            TwoOutputs(), torch.zeros(n_nodes, 1), self.onnx_file)
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        with self.assertRaisesRegex(ValueError, "single output"):
            process.ExecuteFinalizeSolutionStep()


@KratosUnittest.skipUnless(have_torch and have_onnx_export and have_onnxruntime,
                           "Missing required python modules: torch, physicsnemo, onnxruntime.")
class TestOnnxCudaFallbackIsReported(KratosUnittest.TestCase):
    """The silent-fallback guard. Runs on a CPU-only onnxruntime too - that
    is precisely the case it protects."""

    def setUp(self):
        self.onnx_file = Path("test_onnx_cuda_guard.onnx")
        model = torch.nn.Sequential(torch.nn.Linear(3, 3)).double().eval()
        training_utils.ExportOnnxModel(
            model, torch.zeros(1, 3, dtype=torch.float64), str(self.onnx_file))

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.onnx_file))
        KratosUtilities.DeleteFileIfExisting(str(self.onnx_file) + ".card.json")

    def test_RequireDeviceRaisesWhenCudaIsUnavailable(self):
        # A device index that cannot exist: onnxruntime answers with a
        # perfectly working CPU session, which is the whole problem.
        with self.assertRaisesRegex(RuntimeError, "does not exist|onnxruntime-gpu"):
            onnx_utils.CreateOrtSession(self.onnx_file, "cuda:15", require_device=True)

    def test_WithoutRequireDeviceItFallsBackAndStillWorks(self):
        session = onnx_utils.CreateOrtSession(self.onnx_file, "cuda:15")
        self.assertIn("CPUExecutionProvider", session.get_providers())
        # the export fixes the batch dimension at 1
        name = session.get_inputs()[0].name
        values = numpy.zeros((1, 3), dtype=onnx_utils.NumpyDtypeForOrtInput(
            session.get_inputs()[0]))
        self.assertEqual(session.run(None, {name: values})[0].shape[0], 1)


@KratosUnittest.skipUnless(have_torch and have_onnx_export,
                           "Missing required python modules: torch, physicsnemo.")
@KratosUnittest.skipUnless(have_ort_cuda,
                           "Requires an onnxruntime CUDA build (pip install onnxruntime-gpu).")
class TestOnnxGpuInference(KratosUnittest.TestCase):
    """GPU inference. Skips unless onnxruntime itself is a CUDA build."""

    def setUp(self):
        self.onnx_file = Path("test_onnx_gpu.onnx")
        torch.manual_seed(0)
        model = torch.nn.Sequential(
            torch.nn.Linear(4, 8), torch.nn.Tanh(), torch.nn.Linear(8, 2)).eval()
        # the export fixes the batch dimension, so it must match the rows the
        # numerical comparison below feeds
        self.rows = 64
        training_utils.ExportOnnxModel(
            model, torch.zeros(self.rows, 4), str(self.onnx_file))

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.onnx_file))
        KratosUtilities.DeleteFileIfExisting(str(self.onnx_file) + ".card.json")

    def test_SessionActuallyUsesTheCudaProvider(self):
        # get_providers() reports what ORT INSTANTIATED, not what was asked
        # for - the only way to tell a GPU session from a silent CPU one.
        session = onnx_utils.CreateOrtSession(self.onnx_file, "cuda")
        self.assertEqual(session.get_providers()[0], "CUDAExecutionProvider")

    def test_RequireDeviceSucceedsOnAWorkingGpu(self):
        session = onnx_utils.CreateOrtSession(self.onnx_file, "cuda:0", require_device=True)
        self.assertEqual(session.get_providers()[0], "CUDAExecutionProvider")

    def test_GpuMatchesCpuNumerically(self):
        gpu = onnx_utils.CreateOrtSession(self.onnx_file, "cuda")
        cpu = onnx_utils.CreateOrtSession(self.onnx_file, "cpu")
        name = gpu.get_inputs()[0].name
        dtype = onnx_utils.NumpyDtypeForOrtInput(gpu.get_inputs()[0])
        values = numpy.random.default_rng(0).random((self.rows, 4)).astype(dtype)
        on_gpu = gpu.run(None, {name: values})[0]
        on_cpu = cpu.run(None, {name: values})[0]
        # float32 with TF32 matmuls on the GPU: agreement is ~1e-4 relative,
        # not bitwise, so the tolerance is set from the output scale
        scale = max(float(numpy.abs(on_cpu).max()), 1e-12)
        self.assertLess(float(numpy.abs(on_gpu - on_cpu).max()) / scale, 1e-3)


if __name__ == '__main__':
    KratosUnittest.main()
