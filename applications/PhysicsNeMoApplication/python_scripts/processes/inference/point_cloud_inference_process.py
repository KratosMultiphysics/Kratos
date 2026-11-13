"""Process deploying point-cloud models on a model part's nodes.

Point-cloud transformers (physicsnemo.models.transolver.Transolver and
friends) consume per-point features with coordinates, batched as
(1, N, C) - no tessellation, graph or grid required. This process gathers
the nodal input fields plus the (optionally normalized) node coordinates,
runs one forward pass, and writes the (1, N, C_out) prediction back through
the same field-splitting contract as InferenceProcess (whose settings it
extends).

Model interfaces:
- "generic": model(x) with x = (1, N, 3 + C_in), coordinates prepended to
  the features - anything from an MLP to a scripted trunk.
- "transolver": model(fx, embedding) with fx = (1, N, C_in) functional
  features and embedding = (1, N, 3) coordinates, matching
  Transolver.forward.
- "flare": alias of the "transolver" call contract - matching
  physicsnemo.experimental.models.flare.FLARE.forward(fx, embedding)
  (verified signature-compatible). Listed separately for discoverability
  and model-card clarity; FLARE lives in physicsnemo.experimental (no API
  stability guarantee).
- "geotransolver": model(local_embedding, local_positions=..., geometry=...)
  with local_embedding = (1, N, C_in) features, local_positions = (1, N, 3)
  coordinates and geometry = (1, N, 3) coordinates (or None when the model
  was built with geometry_dim=None - set "pass_geometry" to false then),
  matching physicsnemo.experimental.models.geotransolver.GeoTransolver.
  Construct with use_te=False unless transformer_engine is installed.
  Experimental namespace - no API stability guarantee.
- "deeponet": model(x_branch, x_trunk) with x_branch = (1, D) case
  parameters and x_trunk = (N, d) query coordinates, matching
  physicsnemo.experimental.models.xdeeponet.DeepONet in core mode
  (auto_pad=False) with an MLP branch. This is the one interface that is
  NOT pointwise: an operator maps a whole CASE - a conductivity, a load, an
  angle of attack - to a field at arbitrary points, so the per-node input
  fields play no part and the "branch_input" block names the parameters
  instead. It is what RomSurrogateProcess does through a POD basis, without
  the basis. Build the branch as (D_in) -> (width) and the trunk as (d) ->
  (width) with the same width, and set "trunk_dimension" (or leave it at 0
  to follow DOMAIN_SIZE). Experimental namespace - no API stability
  guarantee.
- "globe": model(prediction_points, boundary_meshes, reference_lengths)
  with boundary_meshes a dict of named surface meshes, matching
  physicsnemo.experimental.models.globe.GLOBE. Another non-pointwise
  interface, and a different idea again: the model predicts the interior
  FROM the boundaries, the way an elliptic problem is determined by its
  boundary conditions, so the "globe" block names the sub-model-parts that
  carry the data rather than per-node input fields. Verified to build its
  cluster tree and backpropagate on CPU. Experimental namespace - no API
  stability guarantee.
- "figconvnet": model(vertices, features) with vertices = (1, N, 3) and
  features = (1, N, C_in), matching FIGConvUNet.forward, which returns a
  TUPLE (point features (1, N, C_out), scalar (1, 1) or None). The point
  features are written back as usual; the scalar (a drag-style global
  output) is stashed as ``last_scalar_prediction`` and logged. Notes:
  FIGConvUNet's warp backend is float32-only (the parameter-dtype cast
  covers stock models), construct it with ``has_input_features=True`` and
  ``in_channels`` equal to the total gathered input width, and its default
  aabb of (0,0,0)-(1,1,1) matches ``normalize_coordinates=True``.

Token budgets: a "subsampling" block runs the model on a subset of the
nodes and gives every other node its nearest selected node's prediction -
see utilities.point_subsampling. A transformer's cost grows with the token
count, and a refined Kratos mesh routinely has more nodes than the model
was trained to attend over.

DoMINO remains served through the CAE datapipes (see CaeDatasetExportProcess
and torch_dataset.CreateDoMINODataPipe).

torch is imported lazily on first execution.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.inference_process import InferenceProcess
from KratosMultiphysics.PhysicsNeMoApplication.utilities.nvtx_utils import NvtxRange
from KratosMultiphysics.PhysicsNeMoApplication.utilities import point_subsampling

_MODEL_INTERFACES = ("generic", "transolver", "flare", "geotransolver", "figconvnet",
                     "deeponet", "globe")


def GatherPointCloudCoordinates(model_part, normalize: bool = True,
                                local_only: bool = False):
    """Returns the (N, 3) current node coordinates as a numpy array,
    optionally min-max normalized to [0, 1] per axis (degenerate axes
    left at 0).

    With local_only the owned nodes only (the communicator's LocalMesh) are
    read. Note that normalize then uses each rank's *own* bounding box, so a
    distributed caller wanting rank-independent coordinates must either pass
    normalize=False or normalize against a globally reduced box itself.
    """
    nodes = model_part.GetCommunicator().LocalMesh().Nodes if local_only else model_part.Nodes
    position_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(
        nodes, Kratos.Configuration.Current)
    position_ta.CollectData()
    coordinates = numpy.array(position_ta.data, dtype=numpy.float64)  # (N, 3)
    if normalize:
        low = coordinates.min(axis=0)
        extent = coordinates.max(axis=0) - low
        extent[extent == 0.0] = 1.0  # planar/linear clouds: leave that axis at 0
        coordinates = (coordinates - low) / extent
    return coordinates


def RunPointCloudForward(model, device, model_interface, features, coordinates,
                         pass_geometry: bool = True, enable_grad: bool = False,
                         branch_input=None):
    """One point-cloud forward pass through a model interface.

    Shared by PointCloudInferenceProcess and the CoSimulation surrogate
    wrapper. features/coordinates are (N, C_in)/(N, 3) torch tensors; both
    are cast to the model's parameter dtype before the call. Runs under
    no_grad by default; enable_grad=True keeps the autograd graph (for
    sensitivity computations - see sensitivity_utils).

    branch_input is the (D,) or (B, D) case-parameter vector the "deeponet"
    interface maps from; the other interfaces ignore it.

    Returns:
        (prediction, scalar): the (N, C_out) float64 prediction and the
        figconvnet-style global scalar (None for the other interfaces).
    """
    torch = torch_bridge._TryImportTorch()

    parameter = next(model.parameters(), None)
    dtype = parameter.dtype if parameter is not None else features.dtype
    features = features.to(dtype)
    coordinates = coordinates.to(dtype)
    if branch_input is not None:
        branch_input = branch_input.to(dtype)

    scalar_prediction = None
    grad_context = torch.enable_grad() if enable_grad else torch.no_grad()
    with grad_context, NvtxRange("PhysicsNeMo::Forward"):
        if model_interface in ("transolver", "flare"):
            prediction = model(
                features[None].to(device),
                coordinates[None].to(device))
        elif model_interface == "geotransolver":
            geometry = coordinates[None].to(device) if pass_geometry else None
            prediction = model(
                features[None].to(device),
                local_positions=coordinates[None].to(device),
                geometry=geometry)
        elif model_interface == "figconvnet":
            prediction, scalar = model(
                coordinates[None].to(device),
                features[None].to(device))
            scalar_prediction = (
                float(scalar.reshape(-1)[0]) if scalar is not None else None)
        elif model_interface == "deeponet":
            # an operator, not a pointwise model: the BRANCH takes the case
            # parameters (one vector for the whole solve) and the TRUNK the
            # query coordinates, so the features are not per-point here -
            # branch_input carries them, and every point shares them
            if branch_input is None:
                raise ValueError(
                    "The \"deeponet\" interface needs a branch input - the case "
                    "parameters the operator maps from. Configure the process's "
                    "\"branch_input\" block, or pass branch_input= directly.")
            prediction = model(
                branch_input[None].to(device) if branch_input.ndim == 1
                else branch_input.to(device),
                coordinates.to(device))
        elif model_interface == "generic":  # coordinates prepended to the features
            prediction = model(
                torch.cat([coordinates, features], dim=-1)[None].to(device))
        else:
            raise ValueError(
                f"Unsupported model interface \"{model_interface}\". "
                f"Supported: {', '.join(_MODEL_INTERFACES)}.")
        prediction = prediction.cpu()

    if prediction.ndim != 3 or prediction.shape[0] != 1:
        raise ValueError(
            f"The model must return a (1, N, C_out) prediction; got shape "
            f"{list(prediction.shape)}.")
    return prediction[0].to(torch.float64), scalar_prediction


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "PointCloudInferenceProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return PointCloudInferenceProcess(model, settings["Parameters"])


class PointCloudInferenceProcess(InferenceProcess):
    """Runs point-cloud model inference each output_interval steps."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        # Split the subclass keys off (with their defaults) before the parent
        # validates the shared InferenceProcess settings.
        self.model_interface = "generic"
        self.normalize_coordinates = True
        self.pass_geometry = True  # geotransolver only: forward coordinates as geometry
        self.subsampling_settings = None  # a token budget; see point_subsampling
        # "globe" only: the boundaries whose data drives the interior
        self._globe_settings = None
        if settings.Has("globe"):
            globe = settings["globe"]
            globe.ValidateAndAssignDefaults(Kratos.Parameters("""{
                "boundary_sub_model_parts" : [],
                "boundary_fields"          : [],
                "reference_lengths"        : {},
                "output_names"             : [],
                "source_container"         : "Elements"
            }"""))
            self._globe_settings = globe.Clone()
            settings.RemoveValue("globe")
        # "deeponet" only: the case parameters the operator's branch maps from
        self._branch_process_info_variables = []
        self._branch_properties_variables = []
        self._branch_constants = []
        self._trunk_dimension = 0  # 0 = read DOMAIN_SIZE
        if settings.Has("branch_input"):
            branch = settings["branch_input"]
            branch.ValidateAndAssignDefaults(Kratos.Parameters("""{
                "process_info_variables" : [],
                "properties_variables"   : [],
                "constants"              : []
            }"""))
            self._branch_process_info_variables = [
                branch["process_info_variables"][i].GetString()
                for i in range(branch["process_info_variables"].size())]
            self._branch_properties_variables = [
                branch["properties_variables"][i].GetString()
                for i in range(branch["properties_variables"].size())]
            self._branch_constants = list(branch["constants"].GetVector())
            settings.RemoveValue("branch_input")
        if settings.Has("trunk_dimension"):
            self._trunk_dimension = settings["trunk_dimension"].GetInt()
            settings.RemoveValue("trunk_dimension")
        if settings.Has("subsampling"):
            self.subsampling_settings = settings["subsampling"].Clone()
            settings.RemoveValue("subsampling")
        if settings.Has("model_interface"):
            self.model_interface = settings["model_interface"].GetString()
            settings.RemoveValue("model_interface")
        if settings.Has("normalize_coordinates"):
            self.normalize_coordinates = settings["normalize_coordinates"].GetBool()
            settings.RemoveValue("normalize_coordinates")
        if settings.Has("pass_geometry"):
            self.pass_geometry = settings["pass_geometry"].GetBool()
            settings.RemoveValue("pass_geometry")
        self.last_scalar_prediction = None  # figconvnet-style global output
        self._globe_cache = None
        super().__init__(model, settings)

        if self.model_interface not in _MODEL_INTERFACES:
            raise ValueError(
                f"Unsupported model interface \"{self.model_interface}\". "
                f"Supported: {', '.join(_MODEL_INTERFACES)}.")

    def _GlobeInputs(self):
        """The boundary meshes, reference lengths and output names GLOBE
        takes, rebuilt when the mesh changes and cached otherwise."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import globe_bridge

        if self._globe_settings is None:
            raise ValueError(
                "The \"globe\" interface needs a \"globe\" settings block naming the "
                "boundary sub-model-parts, their fields and the reference lengths.")
        node_count = self.model_part.NumberOfNodes()
        if self._globe_cache is not None and self._globe_cache[0] == node_count:
            return self._globe_cache[1:]

        settings = self._globe_settings
        names = [settings["boundary_sub_model_parts"][i].GetString()
                 for i in range(settings["boundary_sub_model_parts"].size())]
        if not names:
            raise ValueError(
                "\"boundary_sub_model_parts\" is empty: GLOBE predicts the interior "
                "FROM the boundaries, so at least one is required.")
        field_specs, ranks = [], {}
        for i in range(settings["boundary_fields"].size()):
            entry = settings["boundary_fields"][i]
            variable_name = entry["variable_name"].GetString()
            field_specs.append((
                variable_name,
                entry["data_location"].GetString()
                if entry.Has("data_location") else "node_historical"))
            ranks[variable_name] = entry["rank"].GetInt() if entry.Has("rank") else 0
        boundary_meshes = globe_bridge.BuildGlobeBoundaryMeshes(
            self.model_part, names, field_specs,
            source_container=settings["source_container"].GetString(), ranks=ranks)

        reference_lengths = {
            key: settings["reference_lengths"][key].GetDouble()
            for key in settings["reference_lengths"].keys()}
        output_names = [settings["output_names"][i].GetString()
                        for i in range(settings["output_names"].size())]
        if not output_names:
            raise ValueError(
                "\"output_names\" is empty: GLOBE's predictions are named fields, and "
                "the names say which column of the output each one becomes.")
        self._globe_cache = (node_count, boundary_meshes, reference_lengths, output_names)
        return self._globe_cache[1:]

    def _GatherBranchInput(self):
        """The case-parameter vector the "deeponet" branch maps from.

        An operator learns a map from a whole CASE to a field, so its branch
        input is one vector per solve, not one row per node: the
        conductivity and heat flux of a thermal case, the angle of attack of
        an aero case. They are read from wherever Kratos keeps them - the
        ProcessInfo, the elements' Properties, or written in by hand.
        """
        torch = torch_bridge._TryImportTorch()

        values = []
        for name in self._branch_process_info_variables:
            values.append(float(self.model_part.ProcessInfo[
                Kratos.KratosGlobals.GetVariable(name)]))
        if self._branch_properties_variables:
            element = next(iter(self.model_part.Elements), None)
            if element is None:
                raise ValueError(
                    "\"branch_input\" names Properties variables but the model part has "
                    "no elements to read them from.")
            properties = element.Properties
            for name in self._branch_properties_variables:
                values.append(float(properties[Kratos.KratosGlobals.GetVariable(name)]))
        values.extend(self._branch_constants)
        if not values:
            raise ValueError(
                "The \"deeponet\" interface needs a non-empty \"branch_input\" block: "
                "an operator maps case parameters to a field, and no parameters were "
                "named.")
        return torch.tensor(values, dtype=torch.float64)

    def _TrunkCoordinates(self, coordinates):
        """Coordinates cut to the trunk's own dimension."""
        if self._trunk_dimension == 0:
            dimension = int(self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE]) or 3
        else:
            dimension = self._trunk_dimension
        if dimension not in (2, 3):
            raise ValueError(
                f"\"trunk_dimension\" resolved to {dimension}; a DeepONet trunk takes "
                "2-D or 3-D query coordinates.")
        return coordinates[:, :dimension]

    def RunInference(self) -> None:
        torch = torch_bridge._TryImportTorch()

        with NvtxRange("PhysicsNeMo::GatherInputs"):
            features, n_entities = self._GatherFeatures()  # (N, C_in), card-normalized
            coordinates = torch.from_numpy(GatherPointCloudCoordinates(
                self.model_part, self.normalize_coordinates))  # (N, 3)

        # a token budget: the model sees a subset, every node gets a value
        selected = None
        if self.subsampling_settings is not None:
            # selection runs on the UNNORMALIZED coordinates, so a bounding
            # box is written in the model part's own units
            selection_coordinates = GatherPointCloudCoordinates(
                self.model_part, normalize=False)
            selected = point_subsampling.SelectPointSubset(
                selection_coordinates, self.subsampling_settings.Clone())
            if selected is not None:
                index = torch.from_numpy(selected)
                features = features[index]
                coordinates = coordinates[index]
        self._CheckOOD(features)

        branch_input = None
        forward_coordinates = coordinates
        if self.model_interface == "deeponet":
            branch_input = self._GatherBranchInput()
            forward_coordinates = self._TrunkCoordinates(coordinates)

        if self.model_interface == "globe":
            from KratosMultiphysics.PhysicsNeMoApplication.bridges import globe_bridge

            boundary_meshes, reference_lengths, output_names = self._GlobeInputs()
            globe_coordinates = self._TrunkCoordinates(coordinates)

            def forward(model):
                with NvtxRange("PhysicsNeMo::Forward"):
                    return globe_bridge.RunGlobeForward(
                        model, globe_coordinates, boundary_meshes, reference_lengths,
                        output_names)
        else:
            def forward(model):
                prediction, scalar = RunPointCloudForward(
                    model, self._device, self.model_interface, features,
                    forward_coordinates, self.pass_geometry, branch_input=branch_input)
                if self.model_interface == "figconvnet":
                    self.last_scalar_prediction = scalar
                    if scalar is not None:
                        Kratos.Logger.PrintInfo(
                            type(self).__name__,
                            f"figconvnet scalar prediction: {scalar:.6e}")
                return prediction

        prediction, std = self._PredictWithUncertainty(forward)
        if selected is not None:
            prediction = point_subsampling.ExpandToAllPoints(
                prediction, selection_coordinates, selected)
            std = point_subsampling.ExpandToAllPoints(
                std, selection_coordinates, selected)
        with NvtxRange("PhysicsNeMo::WriteOutputs"):
            self._WriteOutputs(prediction, n_entities)
            self._WriteUncertainty(std, n_entities)
