"""Tests for the Triton model-repository export and the remote-inference
client process.

No Triton server runs here, so serving is validated structurally: the
repository layout and generated config.pbtxt are parsed back, the exported
ONNX artifact is round-tripped through ONNX Runtime (including at a DIFFERENT
entity count, which is what proves the dynamic axis took), and the client's
request payload plus write-back are pinned against an injected stub.
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.deployment import triton_export
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import triton_inference_process
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import onnxruntime  # noqa: F401
    have_onnxruntime = have_torch
except ImportError:
    have_onnxruntime = False

try:
    # a tritonclient without its protocol extra raises RuntimeError, not ImportError
    import tritonclient.http  # noqa: F401
    have_tritonclient = True
except (ImportError, RuntimeError):
    have_tritonclient = False


def _ParseConfig(text):
    """Minimal config.pbtxt reader: scalars plus input/output tensor blocks."""
    scalars, tensors = {}, {"input": [], "output": []}
    section, current = None, None
    for raw in text.splitlines():
        line = raw.strip()
        if line in ("input [", "output ["):
            section = line.split()[0]
            continue
        if line == "{":
            current = {}
            continue
        if line == "}":
            if current is not None and section:
                tensors[section].append(current)
            current = None
            continue
        if line == "]":
            section = None
            continue
        if ":" in line:
            key, value = (part.strip() for part in line.split(":", 1))
            value = value.strip('"')
            if value.startswith("[") and value.endswith("]"):
                value = [int(v) for v in value.strip("[]").split(",") if v.strip()]
            target = current if current is not None else scalars
            target[key] = value
    return scalars, tensors


class TestTritonConfigGeneration(KratosUnittest.TestCase):
    """Pure text generation: runs without torch, physicsnemo or Triton."""

    def test_ConfigContents(self):
        text = triton_export.MakeTritonConfig(
            "cavity_surrogate", "onnxruntime_onnx",
            inputs=[{"name": "VELOCITY", "data_type": "TYPE_FP32", "dims": [-1, 3]}],
            outputs=[{"name": "PRESSURE", "data_type": "TYPE_FP32", "dims": [-1, 1]}])
        scalars, tensors = _ParseConfig(text)

        self.assertEqual(scalars["name"], "cavity_surrogate")
        self.assertEqual(scalars["platform"], "onnxruntime_onnx")
        # per-entity tensors have no stackable batch axis
        self.assertEqual(scalars["max_batch_size"], "0")
        self.assertEqual(tensors["input"][0]["name"], "VELOCITY")
        self.assertEqual(tensors["input"][0]["data_type"], "TYPE_FP32")
        self.assertEqual(tensors["input"][0]["dims"], [-1, 3])
        self.assertEqual(tensors["output"][0]["dims"], [-1, 1])

    def test_InstanceGroupAndDynamicBatching(self):
        text = triton_export.MakeTritonConfig(
            "m", "onnxruntime_onnx",
            inputs=[{"name": "in", "data_type": "TYPE_FP32", "dims": [-1, 2]}],
            outputs=[{"name": "out", "data_type": "TYPE_FP32", "dims": [-1, 1]}],
            max_batch_size=8, instance_group=[{"count": 2, "kind": "KIND_GPU"}],
            dynamic_batching=True)
        self.assertIn("instance_group", text)
        self.assertIn("count: 2", text)
        self.assertIn("KIND_GPU", text)
        self.assertIn("dynamic_batching { }", text)

    def test_Validation(self):
        tensor = [{"name": "in", "data_type": "TYPE_FP32", "dims": [-1, 2]}]
        with self.assertRaisesRegex(ValueError, "max_batch_size"):
            triton_export.MakeTritonConfig("m", "p", tensor, tensor, max_batch_size=-1)
        with self.assertRaisesRegex(ValueError, "dynamic_batching"):
            triton_export.MakeTritonConfig("m", "p", tensor, tensor, dynamic_batching=True)
        with self.assertRaisesRegex(ValueError, "at least one input"):
            triton_export.MakeTritonConfig("m", "p", [], tensor)


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestTritonRepositoryExport(KratosUnittest.TestCase):

    def setUp(self):
        self.repository = Path("test_triton_repository")
        torch.manual_seed(0)
        self.model = torch.nn.Sequential(
            torch.nn.Linear(3, 8), torch.nn.Tanh(), torch.nn.Linear(8, 1))
        self.sample = torch.randn(12, 3)

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.repository))
        KratosUtilities.DeleteFileIfExisting("test_triton_card.pt")
        KratosUtilities.DeleteFileIfExisting("test_triton_card.pt.card.json")

    def _Export(self, extra=""):
        return triton_export.ExportTritonModelRepository(
            self.model, self.sample, Kratos.Parameters("""{
                "repository_path" : "%s",
                "model_name"      : "surrogate",
                "model_version"   : 1
                %s
            }""" % (self.repository, extra)))

    def test_RepositoryLayout(self):
        config_file = self._Export()
        self.assertEqual(Path(config_file), self.repository / "surrogate" / "config.pbtxt")
        self.assertTrue((self.repository / "surrogate" / "1" / "model.onnx").is_file())

        scalars, tensors = _ParseConfig(Path(config_file).read_text())
        self.assertEqual(scalars["platform"], "onnxruntime_onnx")
        self.assertEqual(tensors["input"][0]["dims"], [-1, 3])
        self.assertEqual(tensors["output"][0]["dims"], [-1, 1])

    def test_TorchScriptFormat(self):
        config_file = self._Export(', "format": "torchscript"')
        self.assertTrue((self.repository / "surrogate" / "1" / "model.pt").is_file())
        scalars, _ = _ParseConfig(Path(config_file).read_text())
        self.assertEqual(scalars["platform"], "pytorch_libtorch")

    def test_CardDrivesTensorNames(self):
        card_path = "test_triton_card.pt"
        Path(card_path).write_text("stand-in for the checkpoint the card describes")
        model_registry.SaveModelCard(card_path, {
            "input_fields": [{"variable_name": "VELOCITY", "data_location": "node_historical"}],
            "output_fields": [{"variable_name": "PRESSURE", "data_location": "node_historical"}],
        })
        config_file = self._Export(', "card_file": "%s"' % card_path)
        _, tensors = _ParseConfig(Path(config_file).read_text())
        self.assertEqual(tensors["input"][0]["name"], "VELOCITY")
        self.assertEqual(tensors["output"][0]["name"], "PRESSURE")
        # the card travels with the served artifact
        self.assertIsNotNone(model_registry.LoadModelCard(
            str(self.repository / "surrogate" / "1" / "model.onnx")))

    def test_Validation(self):
        with self.assertRaisesRegex(ValueError, "Unsupported format"):
            self._Export(', "format": "tensorrt"')
        with self.assertRaisesRegex(ValueError, "model_version"):
            self._Export(', "model_version": 0')

    @KratosUnittest.skipUnless(have_onnxruntime, "Missing required python module: onnxruntime.")
    def test_ExportedArtifactServesAnyMeshSize(self):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import onnx_utils
        self._Export()
        session = onnx_utils.CreateOrtSession(
            str(self.repository / "surrogate" / "1" / "model.onnx"), "cpu")
        ort_input = session.get_inputs()[0]

        # matches torch on the exported size ...
        with torch.no_grad():
            expected = self.model(self.sample).numpy()
        produced = session.run(None, {ort_input.name: self.sample.numpy()})[0]
        numpy.testing.assert_allclose(produced, expected, rtol=1e-5, atol=1e-6)

        # ... and, thanks to the dynamic entity axis, on a different one too
        other = torch.randn(37, 3)
        with torch.no_grad():
            expected_other = self.model(other).numpy()
        produced_other = session.run(None, {ort_input.name: other.numpy()})[0]
        self.assertEqual(produced_other.shape, (37, 1))
        numpy.testing.assert_allclose(produced_other, expected_other, rtol=1e-5, atol=1e-6)


class _StubResponse:
    def __init__(self, arrays):
        self._arrays = arrays

    def as_numpy(self, name):
        return self._arrays.get(name)


class _StubTritonClient:
    """Records the request and echoes a row-indexed prediction back."""

    def __init__(self, output_name, width=1):
        self.output_name = output_name
        self.width = width
        self.calls = []

    def infer(self, **keywords):
        self.calls.append(keywords)
        rows = keywords["inputs"][0].shape()[0]
        return _StubResponse({self.output_name: numpy.arange(
            rows * self.width, dtype=numpy.float32).reshape(rows, self.width)})


@KratosUnittest.skipUnless(have_torch and have_tritonclient,
                           "Missing required python modules: torch, tritonclient.")
class TestTritonInferenceProcess(KratosUnittest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Remote")
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        for i in range(6):
            node = self.model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.PRESSURE, float(i) + 0.5)

    def _CreateProcess(self, extra=""):
        return triton_inference_process.Factory(Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Remote",
                "model_settings"  : {
                    "url"         : "localhost:8000",
                    "model_name"  : "surrogate",
                    "input_name"  : "VELOCITY",
                    "output_name" : "PRESSURE"
                    %s
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",
                                        "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE",
                                        "data_location" : "node_historical" } ]
            }
        }""" % extra), self.model)

    def test_RequestPayloadAndWriteBack(self):
        process = self._CreateProcess()
        client = _StubTritonClient("PRESSURE")
        process.SetClient(client)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        self.assertEqual(len(client.calls), 1)
        call = client.calls[0]
        self.assertEqual(call["model_name"], "surrogate")

        infer_input = call["inputs"][0]
        self.assertEqual(infer_input.name(), "VELOCITY")
        self.assertEqual(infer_input.shape(), [6, 1])
        self.assertEqual(infer_input.datatype(), "FP32")
        # the gathered field really reached the payload
        sent = numpy.frombuffer(infer_input._get_binary_data(), dtype=numpy.float32)
        numpy.testing.assert_allclose(
            sent, numpy.array([float(i) + 0.5 for i in range(6)], dtype=numpy.float32))
        self.assertEqual(call["outputs"][0].name(), "PRESSURE")

        written = [node.GetSolutionStepValue(Kratos.TEMPERATURE)
                   for node in self.model_part.Nodes]
        self.assertEqual(written, [float(i) for i in range(6)])

    def test_OutputNormalizationFromTheCardIsApplied(self):
        # this process keys its card off "card_file"; the override that
        # makes the lookup work had no test. The stub echoes row indices,
        # so the written field must be std * row + mean.
        card_path = "test_triton_norm_card.pt"
        Path(card_path).write_text("stand-in for the checkpoint the card describes")
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, card_path)
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, card_path + ".card.json")
        mean, std = 30.0, 5.0
        model_registry.SaveModelCard(card_path, {
            "output_normalization": {"type": "mean_std", "mean": [mean], "std": [std]}})

        process = self._CreateProcess(', "card_file": "%s"' % card_path)
        process.SetClient(_StubTritonClient("PRESSURE"))
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        written = [node.GetSolutionStepValue(Kratos.TEMPERATURE)
                   for node in self.model_part.Nodes]
        numpy.testing.assert_allclose(written, [std * i + mean for i in range(6)], rtol=1e-6)

    def test_InputNormalizationFromTheCardReachesThePayload(self):
        # what leaves for the server is the standardized field
        card_path = "test_triton_input_card.pt"
        Path(card_path).write_text("stand-in for the checkpoint the card describes")
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, card_path)
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, card_path + ".card.json")
        mean, std = 2.0, 4.0
        model_registry.SaveModelCard(card_path, {
            "input_normalization": {"type": "mean_std", "mean": [mean], "std": [std]}})

        process = self._CreateProcess(', "card_file": "%s"' % card_path)
        client = _StubTritonClient("PRESSURE")
        process.SetClient(client)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        sent = numpy.frombuffer(client.calls[0]["inputs"][0]._get_binary_data(), dtype=numpy.float32)
        numpy.testing.assert_allclose(
            sent, [((float(i) + 0.5) - mean) / std for i in range(6)], rtol=1e-6)

    def test_ModelVersionAndTimeoutForwarded(self):
        process = self._CreateProcess(', "model_version": "3", "timeout": 2.5')
        client = _StubTritonClient("PRESSURE")
        process.SetClient(client)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        self.assertEqual(client.calls[0]["model_version"], "3")
        self.assertAlmostEqual(client.calls[0]["client_timeout"], 2.5)

    def test_RowMismatchRaises(self):
        process = self._CreateProcess()

        class _WrongRows(_StubTritonClient):
            def infer(self, **keywords):
                return _StubResponse({"PRESSURE": numpy.zeros((3, 1), dtype=numpy.float32)})

        process.SetClient(_WrongRows("PRESSURE"))
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        with self.assertRaisesRegex(ValueError, "expected 6 rows"):
            process.ExecuteFinalizeSolutionStep()

    def test_UnsupportedProtocolRejected(self):
        process = self._CreateProcess(', "protocol": "websocket"')
        process.SetClient(_StubTritonClient("PRESSURE"))
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        with self.assertRaisesRegex(ValueError, "Unsupported protocol"):
            process.ExecuteFinalizeSolutionStep()


if __name__ == '__main__':
    KratosUnittest.main()
