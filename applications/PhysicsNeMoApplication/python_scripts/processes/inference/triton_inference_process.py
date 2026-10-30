"""Process running inference against a Triton Inference Server.

The same gather/split contract as InferenceProcess, but the forward pass is
an RPC to a running Triton server (serving a repository written by
triton_export) instead of a local checkpoint - inference-as-a-service, so
the solver host needs neither torch weights nor a GPU.

"model_settings":
    url: Server address, e.g. "localhost:8000" (http) or "localhost:8001" (grpc).
    model_name: The served model name.
    model_version: "" (default) lets the server pick the latest.
    protocol: "http" (default) or "grpc".
    timeout: Client timeout in seconds (0 = the client's default).
    input_name / output_name: Served tensor names; default to the exporter's
        card-derived names when left empty.
    model_card_policy: "advisory" (default), "strict" or "ignore" - checked
        against "card_file" when given.

tritonclient is an optional runtime dependency, imported lazily. A client can
also be injected with SetClient() - the seam the tests use to pin the request
payload and the write-back without a live server.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.inference_process import InferenceProcess
from KratosMultiphysics.PhysicsNeMoApplication.utilities.nvtx_utils import NvtxRange

_PROTOCOLS = ("http", "grpc")


def _TryImportTritonClient(protocol: str):
    """(client_module, InferInput, InferRequestedOutput) for the protocol."""
    try:
        if protocol == "grpc":
            import tritonclient.grpc as client_module
        else:
            import tritonclient.http as client_module
        return client_module, client_module.InferInput, client_module.InferRequestedOutput
    except (ImportError, RuntimeError) as e:
        # tritonclient installed WITHOUT its protocol extra raises RuntimeError
        # ("the installation does not include http support"), not ImportError
        extra = "grpc" if protocol == "grpc" else "http"
        raise ImportError(
            "PhysicsNeMoApplication.triton_inference_process requires tritonclient with "
            f"its \"{extra}\" extra, which could not be imported. Install it with "
            f"'pip install tritonclient[{extra}]'.") from e


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "TritonInferenceProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return TritonInferenceProcess(model, settings["Parameters"])


class TritonInferenceProcess(InferenceProcess):
    """Runs remote Triton inference each output_interval steps."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__(model, settings)
        # Client creation is deferred to the first execution so that merely
        # constructing the process never requires tritonclient or a server.
        self._client = None
        self._resolved = False

    def SetClient(self, client) -> None:
        """Injects a client (anything exposing Triton's `infer`).

        The seam for testing the request payload and write-back without a
        running server, and for reusing an already-configured client.
        """
        self._client = client

    def _Resolve(self):
        if self._resolved:
            return
        default_settings = Kratos.Parameters("""{
            "url"               : "localhost:8000",
            "model_name"        : "kratos_surrogate",
            "model_version"     : "",
            "protocol"          : "http",
            "timeout"           : 0.0,
            "input_name"        : "input",
            "output_name"       : "output",
            "card_file"         : "",
            "model_card_policy" : "advisory"
        }""")
        self.model_settings.ValidateAndAssignDefaults(default_settings)

        self.protocol = self.model_settings["protocol"].GetString()
        if self.protocol not in _PROTOCOLS:
            raise ValueError(
                f"Unsupported protocol \"{self.protocol}\". Use one of {_PROTOCOLS}.")
        self.model_name = self.model_settings["model_name"].GetString()
        self.model_version = self.model_settings["model_version"].GetString()
        self.input_name = self.model_settings["input_name"].GetString()
        self.output_name = self.model_settings["output_name"].GetString()
        self.timeout = self.model_settings["timeout"].GetDouble()

        policy = self.model_settings["model_card_policy"].GetString()
        card_file = self.model_settings["card_file"].GetString()
        if card_file and policy != "ignore":
            card = model_registry.LoadModelCard(card_file)
            model_registry.ValidateFieldsAgainstCard(
                card, self.input_specs, self.output_specs, type(self).__name__, policy)
        self._resolved = True

    def _GetClient(self):
        self._Resolve()
        if self._client is None:
            client_module, _, _ = _TryImportTritonClient(self.protocol)
            self._client = client_module.InferenceServerClient(
                url=self.model_settings["url"].GetString())
        return self._client

    def _MakeRequest(self, features):
        """Builds the (inputs, outputs) request objects for one feature block."""
        _, InferInput, InferRequestedOutput = _TryImportTritonClient(self.protocol)
        infer_input = InferInput(self.input_name, list(features.shape), "FP32")
        infer_input.set_data_from_numpy(features)
        return [infer_input], [InferRequestedOutput(self.output_name)]

    def _NormalizationCardFile(self):
        """This process keys its model card off "card_file", not
        "checkpoint_file", so the base lookup would find nothing."""
        return self.model_settings["card_file"].GetString() \
            if self.model_settings.Has("card_file") else None

    def RunInference(self) -> None:
        client = self._GetClient()
        torch = torch_bridge._TryImportTorch()

        with NvtxRange("PhysicsNeMo::GatherInputs"):
            features, n_entities = self._GatherFeatures()
        features = numpy.ascontiguousarray(features.numpy().astype(numpy.float32))

        with NvtxRange("PhysicsNeMo::Forward"):
            request_inputs, requested_outputs = self._MakeRequest(features)
            keywords = {"model_name": self.model_name,
                        "inputs": request_inputs,
                        "outputs": requested_outputs}
            if self.model_version:
                keywords["model_version"] = self.model_version
            if self.timeout > 0.0:
                keywords["client_timeout"] = self.timeout
            response = client.infer(**keywords)

        prediction = numpy.asarray(response.as_numpy(self.output_name))
        if prediction is None or prediction.shape[0] != n_entities:
            raise ValueError(
                f"Triton returned {None if prediction is None else prediction.shape} for "
                f"output \"{self.output_name}\"; expected {n_entities} rows. Check that "
                "the served model name and tensor names match the exported repository.")

        with NvtxRange("PhysicsNeMo::WriteOutputs"):
            self._WriteOutputs(
                torch.from_numpy(numpy.ascontiguousarray(prediction.astype(numpy.float64))),
                n_entities)
