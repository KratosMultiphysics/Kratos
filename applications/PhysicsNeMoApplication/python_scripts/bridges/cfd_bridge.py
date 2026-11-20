"""Delegation to physicsnemo-cfd: hybrid initialization and CFD metrics.

Bridges Kratos model parts to ``physicsnemo.cfd`` (the ``nvidia-physicsnemo-
cfd`` package, importable as ``physicsnemo.cfd``):

- **Hybrid initialization**: ``physicsnemo.cfd.hybrid_initialization_tools``
  blends two flow fields (e.g. a potential-flow or ML prediction with a
  reference state) on pyvista meshes. ``CreateFlowfield`` builds a
  ``Flowfield`` from a model part's tessellated surface (or nodal point
  cloud) with the mesh bridge, ``CreateHybridInitialization`` delegates the
  blending, and ``FlowfieldToModelPart`` writes the blended point data back
  onto the nodes.
- **CFD metric registry**: ``physicsnemo.cfd.evaluation.metrics`` registers
  domain-aware aerodynamic metrics (relative L2 of pressure/velocity/shear,
  drag/lift, physics residual norms, UQ calibration metrics, ...).
  ``EvaluateCfdMetrics`` resolves and evaluates them on plain numpy dicts
  keyed by the registry's semantic names ("pressure", "velocity", ...);
  ``ValidationMetricsProcess`` exposes them through its optional
  ``cfd_metrics`` settings block.

physicsnemo-cfd and pyvista are optional runtime dependencies, imported
lazily with actionable error messages.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import GetTensorAdaptor

_NODAL_LOCATIONS = ("node_historical", "node_non_historical")


def _TryImportPhysicsNemoCfd():
    try:
        from physicsnemo.cfd.hybrid_initialization_tools import Flowfield, create_hybrid_initialization
        return Flowfield, create_hybrid_initialization
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.cfd_bridge requires physicsnemo-cfd, which could not be "
            "imported. Install it from source with e.g. 'pip install git+https://github.com/NVIDIA/physicsnemo-cfd' (it is not published to PyPI under any name).") from e


def _TryImportCfdMetrics():
    try:
        from physicsnemo.cfd.evaluation import metrics
        return metrics
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.cfd_bridge requires physicsnemo-cfd, which could not be "
            "imported. Install it from source with e.g. 'pip install git+https://github.com/NVIDIA/physicsnemo-cfd' (it is not published to PyPI under any name).") from e


def _TryImportDominoScaling():
    """physicsnemo-cfd's loader for DoMINO's scaling_factors.pkl.

    The file is pickled referencing a `utils.ScalingFactors` class from the
    DoMINO example repository, which is not an installed module, so a plain
    pickle.load raises ModuleNotFoundError. physicsnemo-cfd ships a
    restricted unpickler that remaps it - never hand-roll a shim for this.
    """
    try:
        # _ScalingUnpickler is private but it is the ONLY restricted loader:
        # the public ScalingFactors.load is a plain pickle.load and raises
        # ModuleNotFoundError on these files, and load_scaling_factors_tensors
        # wraps the unpickler behind a Hydra config object.
        from physicsnemo.cfd.evaluation.models.wrappers.domino.scaling import (
            ScalingFactors, _ScalingUnpickler)
        return ScalingFactors, _ScalingUnpickler
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.cfd_bridge requires physicsnemo-cfd to read DoMINO "
            "scaling factors, which could not be imported. Install it from source with "
            "e.g. 'pip install git+https://github.com/NVIDIA/physicsnemo-cfd' "
            "(it is not published to PyPI under any name).") from e


def _TryImportCfdEvaluationWrappers(name: str = "geotransolver"):
    """Resolves one of physicsnemo-cfd's model evaluation wrapper modules.

    The wrappers (physicsnemo.cfd.evaluation.models.wrappers.*: "geotransolver",
    "geotransolver_gp", "geotransolver_drivaerstar", "domino", "transolver",
    "fignet", "xmgn", ...) are checkpoint/NGC-config driven benchmarking
    recipes; this helper only resolves them so a user can feed them Kratos
    data via ModelPartToPolyData/NodesToPolyData. In-loop deployment goes
    through PointCloudInferenceProcess (model_interface "geotransolver")
    - but note it applies only what the model card says, while these
    wrappers normalize their inputs AND call unscale_model_targets on the
    way out. A pretrained checkpoint deployed in-loop therefore needs its
    scalings expressed in the model card - "output_normalization" for the
    targets, "input_normalization" for the field features - and its
    COORDINATE convention checked separately: upstream centres on the STL
    centre of mass and divides by a reference scale, where
    GatherPointCloudCoordinates min-max normalizes per part instead.
    physicsnemo-cfd is alpha - its wrapper APIs may change.
    """
    try:
        import importlib
        return importlib.import_module(f"physicsnemo.cfd.evaluation.models.wrappers.{name}")
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.cfd_bridge requires physicsnemo-cfd, which could not be "
            "imported. Install it from source with e.g. 'pip install git+https://github.com/NVIDIA/physicsnemo-cfd' (it is not published to PyPI under any name).") from e


def _TryImportPyVista():
    try:
        import pyvista
        return pyvista
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.cfd_bridge requires pyvista, which could not be "
            "imported. Install it with e.g. 'pip install pyvista' (it ships with "
            "nvidia-physicsnemo-cfd).") from e


def _GatherPointData(model_part, provenance, field_specs):
    node_ids = [node.Id for node in model_part.Nodes]
    point_data = {}
    for array_name, variable_name, data_location in field_specs:
        if data_location not in _NODAL_LOCATIONS:
            raise ValueError(
                f"cfd_bridge fields must be nodal ({', '.join(_NODAL_LOCATIONS)}), got "
                f"\"{data_location}\" for \"{variable_name}\".")
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
        data = numpy.array(tensor_adaptor.data, dtype=numpy.float64)
        point_data[array_name] = provenance.GatherNodalField(
            node_ids, data.reshape(len(node_ids), -1))
    return point_data


def ModelPartToPolyData(model_part, field_specs=(), source_container="Conditions",
                        tessellation_mode="smallest_id_diagonal",
                        higher_order_mode="reduce", curved_refinement_levels=2):
    """Tessellates a model part into a pyvista PolyData with point data.

    The triangulation comes from the mesh bridge (watertight smallest-id
    rules; curved mode's synthetic sample points are gathered/interpolated
    like any other point), so the PolyData rows follow the provenance point
    ordering - keep the returned provenance to write results back.

    Args:
        model_part: The model part to tessellate.
        field_specs: iterable of (array_name, variable_name, data_location)
            triples attached as point_data (nodal locations only).
        source_container / tessellation_mode / higher_order_mode /
            curved_refinement_levels: see BuildProvenance.

    Returns:
        (polydata, provenance).
    """
    pyvista = _TryImportPyVista()

    provenance = domain_mesh_builder.BuildProvenance(
        model_part, source_container, tessellation_mode,
        higher_order_mode, curved_refinement_levels)
    if provenance.simplex_cells.shape[1] != 3:
        raise RuntimeError(
            "ModelPartToPolyData needs a surface tessellation (triangles); got simplices "
            f"with {provenance.simplex_cells.shape[1]} nodes - use a surface part or the "
            "Conditions container.")

    n_faces = len(provenance.simplex_cells)
    faces = numpy.empty((n_faces, 4), dtype=numpy.int64)
    faces[:, 0] = 3
    faces[:, 1:] = provenance.simplex_cells
    polydata = pyvista.PolyData(provenance.simplex_points, faces.ravel())
    for name, values in _GatherPointData(model_part, provenance, field_specs).items():
        polydata.point_data[name] = values
    return polydata, provenance


def NodesToPolyData(model_part, field_specs=()):
    """Point-cloud PolyData of the model part's nodes (no faces).

    Args:
        model_part: The model part.
        field_specs: iterable of (array_name, variable_name, data_location).

    Returns:
        (polydata, node_ids): rows follow the Nodes container order.
    """
    pyvista = _TryImportPyVista()

    position_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(
        model_part.Nodes, Kratos.Configuration.Current)
    position_ta.CollectData()
    points = numpy.array(position_ta.data, dtype=numpy.float64).reshape(-1, 3)
    polydata = pyvista.PolyData(points)
    node_ids = [node.Id for node in model_part.Nodes]
    for array_name, variable_name, data_location in field_specs:
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
        polydata.point_data[array_name] = numpy.array(
            tensor_adaptor.data, dtype=numpy.float64).reshape(len(node_ids), -1)
    return polydata, node_ids


_FLOWFIELD_SLOTS = (("velocity_variable", "U"), ("pressure_variable", "p"),
                    ("k_variable", "k"), ("omega_variable", "omega"))


def CreateFlowfield(model_part, settings: Kratos.Parameters):
    """Builds a physicsnemo.cfd Flowfield from a model part.

    Settings (all variables optional - empty strings are skipped):
        {
            "source_container"         : "Conditions",
            "data_location"            : "node_historical",
            "velocity_variable"        : "VELOCITY",
            "pressure_variable"        : "PRESSURE",
            "k_variable"               : "",
            "omega_variable"           : "",
            "fill_missing_with_zeros"  : true,
            "tessellation_mode"        : "smallest_id_diagonal",
            "higher_order_mode"        : "reduce",
            "curved_refinement_levels" : 2
        }

    The hybrid-initialization blend reads all four fields (U, p, k, omega)
    from both flowfields, so with "fill_missing_with_zeros" (default) the
    unset slots become zero arrays; disable it for plain export.

    Returns:
        (flowfield, provenance): the Flowfield uses its default fieldnames
        (U/p/k/omega) as point_data arrays; the provenance aligns the
        PolyData rows for FlowfieldToModelPart.
    """
    Flowfield, _ = _TryImportPhysicsNemoCfd()

    default_settings = Kratos.Parameters("""{
        "source_container"         : "Conditions",
        "data_location"            : "node_historical",
        "velocity_variable"        : "",
        "pressure_variable"        : "",
        "k_variable"               : "",
        "omega_variable"           : "",
        "fill_missing_with_zeros"  : true,
        "tessellation_mode"        : "smallest_id_diagonal",
        "higher_order_mode"        : "reduce",
        "curved_refinement_levels" : 2
    }""")
    settings.ValidateAndAssignDefaults(default_settings)

    data_location = settings["data_location"].GetString()
    field_specs = []
    for setting_name, array_name in _FLOWFIELD_SLOTS:
        variable_name = settings[setting_name].GetString()
        if variable_name:
            field_specs.append((array_name, variable_name, data_location))
    if not field_specs:
        raise ValueError("CreateFlowfield needs at least one flow variable to be set.")

    polydata, provenance = ModelPartToPolyData(
        model_part, field_specs,
        settings["source_container"].GetString(),
        settings["tessellation_mode"].GetString(),
        settings["higher_order_mode"].GetString(),
        settings["curved_refinement_levels"].GetInt())
    if settings["fill_missing_with_zeros"].GetBool():
        for _, array_name in _FLOWFIELD_SLOTS:
            if array_name not in polydata.point_data:
                shape = (polydata.n_points, 3) if array_name == "U" else (polydata.n_points,)
                polydata.point_data[array_name] = numpy.zeros(shape)
    return Flowfield(mesh=polydata), provenance


def _ResolveBlendStrategy(settings, n_rows):
    name = settings["blend_strategy"].GetString()
    if name == "constant":
        weight = settings["constant_weight"].GetDouble()
        if not 0.0 <= weight <= 1.0:
            raise ValueError(f"\"constant_weight\" must be in [0, 1] [ got {weight} ].")
        return lambda flowfield_a, flowfield_b: numpy.full(n_rows, weight)
    if name == "from_field_a_k":
        return None  # upstream default (reads flowfield_a's cell-located k field)
    raise ValueError(
        f"Unsupported \"blend_strategy\" \"{name}\". Supported: \"constant\", \"from_field_a_k\".")


_FIELDNAME_ATTRIBUTES = {"U": "velocity_fieldname", "p": "pressure_fieldname",
                         "k": "k_fieldname", "omega": "omega_fieldname"}


def _BlendOnSharedMesh(Flowfield, flowfield_a, flowfield_b, data_location, blend_fn):
    weight = numpy.asarray(blend_fn(flowfield_a, flowfield_b), dtype=numpy.float64)
    merged = flowfield_a.mesh.copy()
    merged.cell_data.clear()
    merged.point_data.clear()
    merged_data = merged.cell_data if data_location == "cell" else merged.point_data
    for fieldname, attribute in _FIELDNAME_ATTRIBUTES.items():
        blended = None
        for flowfield, w in ((flowfield_a, weight), (flowfield_b, 1.0 - weight)):
            data = flowfield.mesh.cell_data if data_location == "cell" else flowfield.mesh.point_data
            term = numpy.einsum(
                "i...,i->i...",
                numpy.asarray(data[getattr(flowfield, attribute)], dtype=numpy.float64), w)
            blended = term if blended is None else blended + term
        merged_data[fieldname] = blended
    return Flowfield(mesh=merged)


def CreateHybridInitialization(flowfield_a, flowfield_b, settings: Kratos.Parameters = None,
                               blend_strategy=None):
    """Blends two Flowfields via physicsnemo.cfd's hybrid initialization.

    Settings:
        {
            "use_topology_from_mesh"    : "a",
            "flowfield_a_data_location" : "point",
            "flowfield_b_data_location" : "point",
            "blend_strategy"            : "constant",
            "constant_weight"           : 0.5,
            "verbose"                   : false
        }

    "blend_strategy" is "constant" (weight * a + (1 - weight) * b, the only
    strategy that works for point-located data) or "from_field_a_k"
    (upstream's turbulence-threshold default; needs a cell-located k field
    on flowfield_a). A Python callable (flowfield_a, flowfield_b) -> weight
    array can be passed as blend_strategy to override the settings choice.

    When both flowfields live on the same mesh (matching data locations and
    coordinates) the blend is computed here with upstream's exact semantics:
    physicsnemo-cfd 0.0.3a0's no-interpolation branch rebinds both
    flowfields to mesh b and returns flowfield_b unchanged, so the shared
    Kratos-mesh case cannot be delegated. Distinct meshes delegate to
    create_hybrid_initialization (kNN interpolation; needs >= 16 source
    points).

    Returns:
        The blended Flowfield (topology of the chosen input, so the matching
        provenance keeps aligning its points).
    """
    Flowfield, create_hybrid_initialization = _TryImportPhysicsNemoCfd()

    if settings is None:
        settings = Kratos.Parameters("""{}""")
    default_settings = Kratos.Parameters("""{
        "use_topology_from_mesh"     : "a",
        "flowfield_a_data_location"  : "point",
        "flowfield_b_data_location"  : "point",
        "blend_strategy"             : "constant",
        "constant_weight"            : 0.5,
        "verbose"                    : false
    }""")
    settings.ValidateAndAssignDefaults(default_settings)

    location_a = settings["flowfield_a_data_location"].GetString()
    location_b = settings["flowfield_b_data_location"].GetString()

    shared_mesh = (
        location_a == location_b
        and flowfield_a.mesh.n_points == flowfield_b.mesh.n_points
        and numpy.allclose(flowfield_a.mesh.points, flowfield_b.mesh.points))
    n_rows = flowfield_a.mesh.n_cells if location_a == "cell" else flowfield_a.mesh.n_points
    if blend_strategy is None:
        blend_strategy = _ResolveBlendStrategy(settings, n_rows)

    if shared_mesh:
        if blend_strategy is None:  # upstream "from_field_a_k" default
            from physicsnemo.cfd.hybrid_initialization_tools.main import from_field_a_k
            blend_strategy = from_field_a_k
        return _BlendOnSharedMesh(Flowfield, flowfield_a, flowfield_b, location_a, blend_strategy)

    kwargs = {} if blend_strategy is None else {"blend_strategy": blend_strategy}
    return create_hybrid_initialization(
        flowfield_a, flowfield_b,
        use_topology_from_mesh=settings["use_topology_from_mesh"].GetString(),
        flowfield_a_data_location=location_a,
        flowfield_b_data_location=location_b,
        verbose=settings["verbose"].GetBool(),
        **kwargs)


def FlowfieldToModelPart(flowfield, model_part, provenance, field_map,
                         data_location="node_historical") -> None:
    """Writes a Flowfield's point data back onto the model part's nodes.

    The Flowfield's mesh rows must follow the given provenance's point
    ordering (true for CreateFlowfield outputs and for hybrid blends that
    kept that topology). Synthetic (curved) rows are dropped - the write is
    exact on the real nodes, like every scatter in the mesh bridge.

    Args:
        flowfield: The physicsnemo.cfd Flowfield.
        model_part: The model part to write into.
        provenance: The MeshProvenanceMap that built the Flowfield's mesh.
        field_map: {point_data_array_name: kratos_variable_name}.
        data_location: "node_historical" (default) or "node_non_historical".
    """
    for array_name, variable_name in field_map.items():
        values = numpy.asarray(flowfield.mesh.point_data[array_name], dtype=numpy.float64)
        domain_mesh_builder.ScatterFieldBack(
            provenance, values, model_part,
            Kratos.KratosGlobals.GetVariable(variable_name), data_location)


def ListCfdMetrics():
    """Names of the metrics registered by physicsnemo-cfd."""
    return sorted(_TryImportCfdMetrics().list_metrics())


def EvaluateCfdMetrics(metric_specs, ground_truth: dict, predictions: dict) -> dict:
    """Evaluates registered physicsnemo-cfd metrics on numpy field dicts.

    Args:
        metric_specs: iterable of (name, domain) pairs - domain "surface"
            or "volume" per the registry.
        ground_truth / predictions: dicts of numpy arrays keyed by the
            registry's SEMANTIC names ("pressure", "velocity", ...) - each
            metric documents which keys it reads; missing keys yield NaN in
            the registry's numpy fallback path.

    Returns:
        {name: float}.
    """
    metrics = _TryImportCfdMetrics()
    values = {}
    for name, domain in metric_specs:
        metric = metrics.get_metric(name, domain=domain)
        values[name] = float(metric(ground_truth, predictions))
    return values
