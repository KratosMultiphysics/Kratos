"""Process running an exported ONNX model inside a Kratos solution loop.

Same gather/split contract as InferenceProcess (whose settings it extends),
but the forward pass runs through a cached onnxruntime InferenceSession
instead of torch - the deployment artifact is a portable .onnx file
(produced by training_utils.ExportOnnxModel) and physicsnemo is not needed
at inference time. The model card sidecar ("<onnx_file>.card.json") is
validated exactly like a checkpoint's.

"model_settings":
    onnx_file: Path to the .onnx model.
    device: "cpu" (default) or "cuda"[:index] - ORT execution provider.
    model_card_policy: "advisory" (default), "strict" or "ignore".

The exported graph must take one (n_entities, total_input_width) tensor and
return one (n_entities, total_output_width) tensor; inputs are cast to the
graph's input dtype (float32 for models exported from float32 weights).

torch is still used as the array bridge to Kratos (lazily, on first
execution); onnxruntime replaces it only for the forward pass.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.deployment import onnx_utils
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.inference_process import InferenceProcess
from KratosMultiphysics.PhysicsNeMoApplication.utilities.nvtx_utils import NvtxRange


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "OnnxInferenceProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return OnnxInferenceProcess(model, settings["Parameters"])


class OnnxInferenceProcess(InferenceProcess):
    """Runs ONNX model inference each output_interval steps."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__(model, settings)
        # Session creation is deferred to the first execution so that merely
        # constructing the process never requires onnxruntime.
        self._session = None

    def _GetSession(self):
        if self._session is None:
            default_settings = Kratos.Parameters("""{
                "onnx_file"         : "PLEASE_SPECIFY_ONNX_FILE",
                "device"            : "cpu",
                "require_device"    : false,
                "model_card_policy" : "advisory"
            }""")
            self.model_settings.ValidateAndAssignDefaults(default_settings)

            onnx_file = self.model_settings["onnx_file"].GetString()
            policy = self.model_settings["model_card_policy"].GetString()
            if policy != "ignore":
                card = model_registry.LoadModelCard(onnx_file)
                model_registry.ValidateFieldsAgainstCard(
                    card, self.input_specs, self.output_specs, type(self).__name__, policy)
            self._session = onnx_utils.CreateOrtSession(
                onnx_file, self.model_settings["device"].GetString(),
                require_device=self.model_settings["require_device"].GetBool())
        return self._session

    def _NormalizationCardFile(self):
        """This process keys its model card off "onnx_file", not
        "checkpoint_file", so the base lookup would find nothing."""
        return self.model_settings["onnx_file"].GetString() \
            if self.model_settings.Has("onnx_file") else None

    def RunInference(self) -> None:
        session = self._GetSession()
        torch = torch_bridge._TryImportTorch()

        with NvtxRange("PhysicsNeMo::GatherInputs"):
            features, n_entities = self._GatherFeatures()

        ort_inputs = session.get_inputs()
        if len(ort_inputs) != 1:
            raise ValueError(
                f"The ONNX model must take a single input tensor; got {len(ort_inputs)} "
                f"({[i.name for i in ort_inputs]}). Export a model whose forward takes one "
                "(n_entities, total_width) tensor.")
        ort_input = ort_inputs[0]
        dtype = onnx_utils.NumpyDtypeForOrtInput(ort_input)
        features = features.numpy().astype(dtype)

        with NvtxRange("PhysicsNeMo::Forward"):
            outputs = session.run(None, {ort_input.name: features})
        if len(outputs) != 1:
            raise ValueError(
                f"The ONNX model must return a single output tensor; got {len(outputs)}. "
                "Export a model whose forward returns one (n_entities, total_width) tensor.")
        prediction = torch.from_numpy(
            numpy.ascontiguousarray(outputs[0]).astype(numpy.float64))

        with NvtxRange("PhysicsNeMo::WriteOutputs"):
            self._WriteOutputs(prediction, n_entities)
