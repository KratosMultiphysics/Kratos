"""Process running a graph neural network on the mesh graph each step.

Deploys MeshGraphNet-style models inside a Kratos solution loop. The graph
topology is extracted once (ExecuteInitialize); node features are re-gathered
per execution; edge features are recomputed per execution only when
"update_edge_features" is enabled (moving meshes).

Model interfaces ("model_interface" setting):

- "meshgraphnet" (default): forward(node_features, edge_features, graph)
  - physicsnemo MeshGraphNet and anything sharing its signature.
- "meshgraphkan": the same call; MeshGraphKAN swaps the node encoder for a
  Fourier Kolmogorov-Arnold layer (extra ctor arg num_harmonics; note it
  ignores num_layers_node_encoder/hidden_dim_node_encoder). Named separately
  only so a deployment reads as what it is.
- "bistride": forward(node_features, edge_features, graph, ms_edges, ms_ids)
  - BiStrideMeshGraphNet runs a U-Net over coarsened copies of the mesh
  graph. It requires graph.pos (no fallback) and the multiscale tables built
  by graph_bridge.BuildBistrideHierarchy; "multiscale_levels" must match the
  model's num_mesh_levels, and bistride_pos_dim must match the coordinate
  width (3 here).
- "hybrid": forward(node_features, mesh_edge_features, world_edge_features,
  graph) - HybridMeshGraphNet adds proximity "world" edges (the "world_edges"
  connectivity block) to the topological mesh edges. The two edge sets live
  in ONE graph object, concatenated mesh-first, and the model splits them
  back positionally by row count; both feature sets must share the model's
  single input_dim_edges.

torch/torch_geometric/physicsnemo are imported lazily.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.deployment import ood_guard_utils
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.utilities.nvtx_utils import NvtxRange

_EXECUTION_POINTS = ("initialize_solution_step", "finalize_solution_step")
_MODEL_INTERFACES = ("meshgraphnet", "meshgraphkan", "bistride", "hybrid")


def RunGraphForward(model, device, model_interface: str, node_features, edge_features,
                    graph, multiscale=None, world_edge_features=None):
    """Runs one graph model forward, dispatching on the model interface.

    Shared by GraphInferenceProcess and any caller holding an already-built
    graph (the same split point as point_cloud_inference_process's
    RunPointCloudForward).

    Args:
        model: The loaded torch model.
        device: Target torch device.
        model_interface: One of _MODEL_INTERFACES.
        node_features: (N, F) array.
        edge_features: (E_mesh, 4) array - MESH edges only for "hybrid".
        graph: The torch_geometric graph (already on device; carrying pos for
            "bistride" and the concatenated mesh+world edges for "hybrid").
        multiscale: (ms_edges, ms_ids) device tensors, "bistride" only.
        world_edge_features: (E_world, 4) array, "hybrid" only.

    Returns:
        (N, C_out) float64 numpy prediction.
    """
    torch = torch_bridge._TryImportTorch()

    parameter = next(model.parameters(), None)
    dtype = parameter.dtype if parameter is not None else torch.float64
    node_tensor = torch.from_numpy(node_features).to(device, dtype)
    edge_tensor = torch.from_numpy(edge_features).to(device, dtype)

    with torch.no_grad(), NvtxRange("PhysicsNeMo::Forward"):
        if model_interface in ("meshgraphnet", "meshgraphkan"):
            prediction = model(node_tensor, edge_tensor, graph)
        elif model_interface == "bistride":
            if multiscale is None:
                raise ValueError(
                    "The \"bistride\" interface needs the multiscale tables; build them "
                    "with graph_bridge.BuildBistrideHierarchy.")
            ms_edges, ms_ids = multiscale
            prediction = model(node_tensor, edge_tensor, graph, ms_edges, ms_ids)
        elif model_interface == "hybrid":
            if world_edge_features is None:
                raise ValueError(
                    "The \"hybrid\" interface needs world edge features; build them with "
                    "graph_bridge.BuildWorldEdges.")
            prediction = model(node_tensor, edge_tensor,
                               torch.from_numpy(world_edge_features).to(device, dtype), graph)
        else:
            raise ValueError(
                f"Unsupported model interface \"{model_interface}\". "
                f"Supported: {', '.join(_MODEL_INTERFACES)}.")
    return prediction.cpu().double().numpy()


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "GraphInferenceProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return GraphInferenceProcess(model, settings["Parameters"])


class GraphInferenceProcess(Kratos.Process):
    """Runs GNN inference on the mesh graph every output_interval steps."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        # optional OOD guard on the gathered node features (see ood_guard_utils)
        ood_settings = Kratos.Parameters("{}")
        if settings.Has("ood_guard"):
            ood_settings = settings["ood_guard"].Clone()
            settings.RemoveValue("ood_guard")
        self._ood_guard = ood_guard_utils.GuardCheck(ood_settings)

        default_settings = Kratos.Parameters("""{
            "model_part_name"      : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "model_settings"       : {},
            "input_fields"         : [
                {
                    "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "data_location" : "node_historical"
                }
            ],
            "output_fields"        : [
                {
                    "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "data_location" : "node_historical"
                }
            ],
            "model_interface"      : "meshgraphnet",
            "source_container"     : "Elements",
            "update_edge_features" : false,
            "multiscale_levels"    : 1,
            "world_edges"          : {
                "type"          : "radius",
                "radius"        : 0.1,
                "max_neighbors" : 16,
                "backend"       : "auto"
            },
            "execution_point"      : "finalize_solution_step",
            "output_interval"      : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        for key in ("input_fields", "output_fields"):
            for i in range(settings[key].size()):
                settings[key][i].ValidateAndAssignDefaults(default_settings[key][0])

        self.model_part = model[settings["model_part_name"].GetString()]
        self.model_settings = settings["model_settings"].Clone()
        self.input_specs = [
            (settings["input_fields"][i]["variable_name"].GetString(),
             settings["input_fields"][i]["data_location"].GetString())
            for i in range(settings["input_fields"].size())
        ]
        self.output_specs = [
            (settings["output_fields"][i]["variable_name"].GetString(),
             settings["output_fields"][i]["data_location"].GetString())
            for i in range(settings["output_fields"].size())
        ]
        self.model_interface = settings["model_interface"].GetString()
        if self.model_interface not in _MODEL_INTERFACES:
            raise ValueError(
                f"Unsupported model interface \"{self.model_interface}\". "
                f"Supported: {', '.join(_MODEL_INTERFACES)}.")
        self.source_container = settings["source_container"].GetString()
        self.update_edge_features = settings["update_edge_features"].GetBool()
        self.multiscale_levels = settings["multiscale_levels"].GetInt()
        if self.multiscale_levels < 1:
            raise ValueError(
                f"\"multiscale_levels\" must be >= 1 [ multiscale_levels = {self.multiscale_levels} ].")
        self.world_edges = settings["world_edges"].Clone()
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
        self._graph = None
        self._edge_index = None
        self._edge_features = None
        self._node_ids = None
        self._multiscale = None
        self._world_edge_features = None
        self._scatter_rows = None
        self._normalization = None
        self._input_normalization = None

    def ExecuteInitialize(self) -> None:
        # Topology is fixed; extract the graph once.
        _, self._edge_index, self._edge_features, self._node_ids = graph_bridge.BuildGraph(
            self.model_part, (), self.source_container)

    def ExecuteInitializeSolutionStep(self) -> None:
        if self.execution_point == "initialize_solution_step":
            self._ExecuteIfDue()

    def ExecuteFinalizeSolutionStep(self) -> None:
        if self.execution_point == "finalize_solution_step":
            self._ExecuteIfDue()

    def _ExecuteIfDue(self) -> None:
        if self.model_part.ProcessInfo[Kratos.STEP] % self.output_interval == 0:
            self.RunInference()

    def RunInference(self) -> None:
        torch = torch_bridge._TryImportTorch()
        if self._edge_index is None:
            self.ExecuteInitialize()
        elif len(self._node_ids) != self.model_part.NumberOfNodes():
            # topology is extracted once, but AdaptiveRemeshProcess ships in
            # this application, so the mesh can change under a running
            # process; a stale edge index silently describes the old mesh
            self._graph = None
            self._scatter_rows = None
            self.ExecuteInitialize()
        if self._model is None:
            self._model, self._device = model_registry.LoadModelWithCardCheck(
                self.model_settings, self.input_specs, self.output_specs, type(self).__name__)
            self._normalization = model_registry.LoadOutputNormalization(self.model_settings)
            self._input_normalization = model_registry.LoadInputNormalization(self.model_settings)
        if self._graph is None:
            self._BuildGraphObject()

        # only the node features change per step; the edge index is topology
        # and was extracted once in ExecuteInitialize. Re-running BuildGraph
        # here would rebuild and discard the entire edge set - three orders of
        # magnitude more than the gather it is wanted for.
        with NvtxRange("PhysicsNeMo::GatherNodeFeatures"):
            node_features = graph_bridge.GatherNodeFeatures(
                self.model_part, self.input_specs, len(self._node_ids))
            # the card's "input_normalization": what the model saw in training
            node_features = model_registry.ApplyInputNormalization(
                node_features, self._input_normalization)
        if self.update_edge_features:
            # geometry, not topology: recomputed from the cached index
            with NvtxRange("PhysicsNeMo::ComputeEdgeFeatures"):
                edge_features = graph_bridge.ComputeEdgeFeatures(
                    self.model_part, self._edge_index)
        else:
            edge_features = self._edge_features
        if self._ood_guard.enabled:
            self._ood_guard.Check(torch.from_numpy(node_features), type(self).__name__)

        prediction = RunGraphForward(
            self._model, self._device, self.model_interface, node_features, edge_features,
            self._graph, self._multiscale, self._world_edge_features)
        # a model trained on normalized targets predicts normalized values;
        # the card's "output_normalization" makes them physical
        prediction = model_registry.ApplyOutputNormalization(prediction, self._normalization)

        with NvtxRange("PhysicsNeMo::ScatterNodeFeatures"):
            if self._scatter_rows is None:
                self._scatter_rows = graph_bridge.BuildScatterRows(
                    self.model_part, self._node_ids)
            graph_bridge.ScatterNodeFeatures(
                self.model_part, self._node_ids, prediction, self.output_specs,
                rows=self._scatter_rows)

    def _BuildGraphObject(self) -> None:
        """Builds the (cached, device-resident) graph object this interface needs.

        Everything the model reads besides the per-step features lives here:
        node coordinates for bistride, the multiscale tables, and the world
        edges appended after the mesh edges for hybrid. All of it is moved to
        the model's device once - including ms_ids, which the model itself
        does NOT move (it only moves ms_edges).
        """
        torch = torch_bridge._TryImportTorch()

        parameter = next(self._model.parameters(), None)
        dtype = parameter.dtype if parameter is not None else torch.float64

        positions = None
        edge_index = self._edge_index
        if self.model_interface == "bistride":
            positions = graph_bridge.NodePositions(self.model_part)
        elif self.model_interface == "hybrid":
            world_edge_index, self._world_edge_features = graph_bridge.BuildWorldEdges(
                self.model_part, self.world_edges)
            edge_index = graph_bridge.ConcatenateEdgeSets(self._edge_index, world_edge_index)

        self._graph = graph_bridge.ToPyGGraph(
            edge_index, len(self._node_ids), positions).to(self._device)
        if self._graph.pos is not None:
            # pos feeds the bistride edge MLPs, so it must match their dtype
            self._graph.pos = self._graph.pos.to(dtype)

        if self.model_interface == "bistride":
            ms_edges, ms_ids = graph_bridge.BuildBistrideHierarchy(
                self._edge_index, len(self._node_ids),
                graph_bridge.NodePositions(self.model_part), self.multiscale_levels)
            self._multiscale = (
                [torch.from_numpy(edges).to(self._device) for edges in ms_edges],
                [torch.from_numpy(ids).to(self._device) for ids in ms_ids])
