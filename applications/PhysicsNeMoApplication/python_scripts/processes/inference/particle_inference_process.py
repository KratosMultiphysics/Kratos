"""Process deploying Learning-to-Simulate particle surrogates.

The Lagrangian counterpart of GraphInferenceProcess: each due step, the
proximity graph of the model part's nodes is REBUILT from the current
positions (particle_bridge - connectivity is proximity, not elements), the
standard Learning-to-Simulate features are gathered (velocity history,
oldest first, plus an optional node-type one-hot), and the model predicts
per-particle ACCELERATION. The state is advanced with semi-implicit Euler
using ProcessInfo[DELTA_TIME]:

    v_new = v + dt * a        (written into VELOCITY)
    x_new = x + dt * v_new    (nodes are MOVED; DISPLACEMENT = x_new - X0)

The velocity history is kept in the process itself (a rolling deque of the
gathered VELOCITY states, like TimeSeriesInferenceProcess - NOT the
historical buffer, whose freshly-cloned slot duplicates the current value
in a self-driving loop): with history_size K > 1, the first K-1 due steps
only warm the history up (logged, nothing advanced). The windows then match
CreateParticleTrajectoryDataset's training layout exactly.

Model interfaces:
- "meshgraphnet" (default): model(node_features, edge_features, graph) with
  a torch_geometric graph (graph_bridge.ToPyGGraph) - the
  physicsnemo.models.meshgraphnet.MeshGraphNet contract (construct with
  input_dim_edges=4, output_dim=3).
- "tensor": model(node_features, edge_features, edge_index) with plain
  tensors - scriptable stubs and custom GNNs, no torch_geometric needed.

torch is imported lazily on first execution.
"""

import collections

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.deployment import ood_guard_utils
from KratosMultiphysics.PhysicsNeMoApplication.bridges import particle_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.utilities.nvtx_utils import NvtxRange

_EXECUTION_POINTS = ("initialize_solution_step", "finalize_solution_step")
_MODEL_INTERFACES = ("meshgraphnet", "tensor")


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "ParticleInferenceProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return ParticleInferenceProcess(model, settings["Parameters"])


class ParticleInferenceProcess(Kratos.Process):
    """Autoregressive particle-dynamics surrogate on a model part's nodes."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        # optional OOD guard on the node features (see ood_guard_utils)
        ood_settings = Kratos.Parameters("{}")
        if settings.Has("ood_guard"):
            ood_settings = settings["ood_guard"].Clone()
            settings.RemoveValue("ood_guard")
        self._ood_guard = ood_guard_utils.GuardCheck(ood_settings)

        default_settings = Kratos.Parameters("""{
            "model_part_name"    : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "model_settings"     : {},
            "model_interface"    : "meshgraphnet",
            "connectivity"       : {},
            "history_size"       : 2,
            "node_type_variable" : "",
            "num_node_types"     : 0,
            "execution_point"    : "finalize_solution_step",
            "output_interval"    : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.model_settings = settings["model_settings"].Clone()
        self.connectivity = settings["connectivity"].Clone()
        self.model_interface = settings["model_interface"].GetString()
        if self.model_interface not in _MODEL_INTERFACES:
            raise ValueError(
                f"Unsupported model interface \"{self.model_interface}\". "
                f"Supported: {', '.join(_MODEL_INTERFACES)}.")
        self.history_size = settings["history_size"].GetInt()
        if self.history_size < 1:
            raise ValueError(f"\"history_size\" must be >= 1 [ history_size = {self.history_size} ].")
        self.node_type_variable = settings["node_type_variable"].GetString()
        self.num_node_types = settings["num_node_types"].GetInt()
        if self.node_type_variable and self.num_node_types < 2:
            raise ValueError(
                "\"num_node_types\" must be >= 2 when \"node_type_variable\" is set "
                f"[ num_node_types = {self.num_node_types} ].")

        self.execution_point = settings["execution_point"].GetString()
        if self.execution_point not in _EXECUTION_POINTS:
            raise ValueError(
                f"Unsupported execution point \"{self.execution_point}\". "
                f"Supported: {', '.join(_EXECUTION_POINTS)}.")
        self.output_interval = settings["output_interval"].GetInt()
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")

        self._model = None
        self._device = None
        self._history = collections.deque(maxlen=self.history_size)

    def ExecuteInitializeSolutionStep(self) -> None:
        if self.execution_point == "initialize_solution_step":
            self._ExecuteIfDue()

    def ExecuteFinalizeSolutionStep(self) -> None:
        if self.execution_point == "finalize_solution_step":
            self._ExecuteIfDue()

    def _ExecuteIfDue(self) -> None:
        if self.model_part.ProcessInfo[Kratos.STEP] % self.output_interval != 0:
            return
        self._history.append(self._GatherVelocities())
        if len(self._history) < self.history_size:
            Kratos.Logger.PrintInfo(
                type(self).__name__,
                f"Warming up the velocity history ({len(self._history)}/"
                f"{self.history_size}); no prediction yet.")
            return
        self.RunInference()

    def _GatherVelocities(self):
        velocities = numpy.empty((self.model_part.NumberOfNodes(), 3))
        for row, node in enumerate(self.model_part.Nodes):
            velocities[row] = node.GetSolutionStepValue(Kratos.VELOCITY)
        return velocities

    def _GatherNodeFeatures(self):
        features = numpy.concatenate(list(self._history), axis=1)  # oldest first
        # The card's "input_normalization" standardizes the velocity history
        # exactly as CreateParticleTrajectoryDataset(normalize=True) did:
        # width K*3, BEFORE the node-type one-hot is appended. Fed raw, a
        # model trained that way drifted 18% over a rollout.
        features = model_registry.ApplyInputNormalization(
            features, getattr(self, "_input_normalization", None))
        if self.node_type_variable:
            variable = Kratos.KratosGlobals.GetVariable(self.node_type_variable)
            one_hot = numpy.zeros((features.shape[0], self.num_node_types))
            for row, node in enumerate(self.model_part.Nodes):
                node_type = int(node.GetValue(variable))
                if not 0 <= node_type < self.num_node_types:
                    raise ValueError(
                        f"Node {node.Id} has type {node_type}, outside "
                        f"[0, {self.num_node_types}).")
                one_hot[row, node_type] = 1.0
            features = numpy.concatenate([features, one_hot], axis=1)
        return features

    def RunInference(self) -> None:
        torch = torch_bridge._TryImportTorch()
        if self._model is None:
            input_specs = [("VELOCITY", "node_historical")]
            self._model, self._device = model_registry.LoadModelWithCardCheck(
                self.model_settings, input_specs,
                [("ACCELERATION", "node_historical")], type(self).__name__)
            # A model trained on standardized accelerations predicts
            # standardized accelerations. This output is integrated TWICE
            # (v += dt*a, then x += dt*v), so leaving it normalized
            # compounds the error straight into node positions.
            self._normalization = model_registry.LoadOutputNormalization(
                self.model_settings)
            self._input_normalization = model_registry.LoadInputNormalization(
                self.model_settings)

        with NvtxRange("PhysicsNeMo::BuildParticleGraph"):
            node_features, edge_index, edge_features, node_ids = particle_bridge.BuildParticleGraph(
                self.model_part, self.connectivity)
            node_features = self._GatherNodeFeatures()
        if self._ood_guard.enabled:
            self._ood_guard.Check(torch.from_numpy(node_features), type(self).__name__)

        parameter = next(self._model.parameters(), None)
        dtype = parameter.dtype if parameter is not None else torch.float64
        with torch.no_grad(), NvtxRange("PhysicsNeMo::Forward"):
            if self.model_interface == "meshgraphnet":
                graph = graph_bridge.ToPyGGraph(edge_index, len(node_ids)).to(self._device)
                acceleration = self._model(
                    torch.from_numpy(node_features).to(self._device, dtype),
                    torch.from_numpy(edge_features).to(self._device, dtype),
                    graph).cpu().double().numpy()
            else:  # tensor: plain (nodes, edges, edge_index) tensors
                acceleration = self._model(
                    torch.from_numpy(node_features).to(self._device, dtype),
                    torch.from_numpy(edge_features).to(self._device, dtype),
                    torch.from_numpy(edge_index).to(self._device)).cpu().double().numpy()

        acceleration = model_registry.ApplyOutputNormalization(
            acceleration, getattr(self, "_normalization", None))
        if tuple(acceleration.shape) != (len(node_ids), 3):
            raise ValueError(
                f"The model must return an ({len(node_ids)}, 3) acceleration; got shape "
                f"{acceleration.shape}.")

        with NvtxRange("PhysicsNeMo::AdvanceParticles"):
            self._AdvanceState(acceleration)

    def _AdvanceState(self, acceleration) -> None:
        """Semi-implicit Euler update; writes ACCELERATION/VELOCITY/DISPLACEMENT
        and moves the nodes."""
        delta_time = self.model_part.ProcessInfo[Kratos.DELTA_TIME]
        if delta_time <= 0.0:
            raise ValueError(
                f"ProcessInfo[DELTA_TIME] must be > 0 to advance particles "
                f"[ delta_time = {delta_time} ].")
        for row, node in enumerate(self.model_part.Nodes):
            a = acceleration[row]
            v = numpy.array(node.GetSolutionStepValue(Kratos.VELOCITY)) + delta_time * a
            node.SetSolutionStepValue(Kratos.ACCELERATION, a.tolist())
            node.SetSolutionStepValue(Kratos.VELOCITY, v.tolist())
            node.X += delta_time * v[0]
            node.Y += delta_time * v[1]
            node.Z += delta_time * v[2]
            node.SetSolutionStepValue(
                Kratos.DISPLACEMENT, [node.X - node.X0, node.Y - node.Y0, node.Z - node.Z0])
