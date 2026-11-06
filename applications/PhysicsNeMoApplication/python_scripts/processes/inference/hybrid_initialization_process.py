"""Process warm-starting a solve from a trained model's prediction.

Runs one inference in ExecuteBeforeSolutionLoop and writes the predicted
fields (typically into historical solution-step variables) so the first
nonlinear solve starts from the ML prediction instead of a cold state,
reducing iteration counts — the "hybrid initialization" idea popularized by
NVIDIA's physicsnemo-cfd. This implementation performs a plain forward pass
through the model loaded via model_registry; to blend the prediction with a
second flow state through physicsnemo-cfd's own recipes, use cfd_bridge
(CreateFlowfield / CreateHybridInitialization / FlowfieldToModelPart).

Reuses all the gather/forward/scatter machinery of InferenceProcess.
"""

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.inference_process import InferenceProcess


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "HybridInitializationProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return HybridInitializationProcess(model, settings["Parameters"])


class HybridInitializationProcess(InferenceProcess):
    """Runs inference exactly once, before the solution loop starts."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        # The inference machinery is reused with its execution points
        # disabled: this process only ever runs in ExecuteBeforeSolutionLoop.
        settings = settings.Clone()
        if settings.Has("execution_point"):
            raise ValueError(
                "\"execution_point\" is not a valid setting of HybridInitializationProcess: "
                "it always executes once, before the solution loop.")
        if settings.Has("output_interval"):
            raise ValueError(
                "\"output_interval\" is not a valid setting of HybridInitializationProcess: "
                "it always executes exactly once.")
        settings.AddString("execution_point", "finalize_solution_step")
        settings.AddInt("output_interval", 1)
        super().__init__(model, settings)
        self._initialized = False

    def ExecuteBeforeSolutionLoop(self) -> None:
        if not self._initialized:
            self.RunInference()
            self._initialized = True

    def ExecuteInitializeSolutionStep(self) -> None:
        pass  # deliberately inert: this process only warm-starts

    def ExecuteFinalizeSolutionStep(self) -> None:
        pass  # deliberately inert: this process only warm-starts
