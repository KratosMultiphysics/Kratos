"""Alignment between physicsnemo's DistributedManager and Kratos MPI.

physicsnemo.distributed.DistributedManager and Kratos's DataCommunicator are
two independent views of the same set of processes. This module initializes
the torch/physicsnemo side consistently with the Kratos side and fails loudly
when the two disagree (e.g. torch picked up a different launcher's
environment variables than the MPI ranks Kratos runs under).

Scope note: this covers the "one job, both worlds" case — e.g. distributed
inference inside an MPI-parallel analysis. For active learning, keeping
training and solving in separate OS processes (the SubprocessBackend) remains
the recommended architecture.

torch/physicsnemo are imported lazily.
"""

import os

import KratosMultiphysics as Kratos

_LAUNCHER_ENV_VARIABLES = ("RANK", "SLURM_PROCID", "OMPI_COMM_WORLD_RANK")


def _TryImportDistributedManager():
    try:
        from physicsnemo.distributed.manager import DistributedManager
        return DistributedManager
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.distributed_utils requires physicsnemo, which could not "
            "be imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e


def _CheckConsistency(manager, data_communicator) -> None:
    """Raises when torch's and Kratos's views of the job disagree."""
    if manager.rank != data_communicator.Rank() or manager.world_size != data_communicator.Size():
        raise RuntimeError(
            "physicsnemo's DistributedManager and Kratos's DataCommunicator disagree: "
            f"torch sees rank {manager.rank} of {manager.world_size}, Kratos sees rank "
            f"{data_communicator.Rank()} of {data_communicator.Size()}. This usually means "
            "torch.distributed was initialized from a different launcher environment than "
            "the MPI ranks Kratos runs under.")


def InitializeDistributedManager(data_communicator=None,
                                 addr: str = "localhost",
                                 port: str = "12355",
                                 backend=None):
    """Initializes physicsnemo's DistributedManager consistently with Kratos.

    - Serial Kratos run: explicit single-process setup.
    - Distributed run with a recognized launcher environment (torch env vars,
      SLURM, OpenMPI): physicsnemo's own auto-detection, then a consistency
      check against the DataCommunicator.
    - Distributed run without such an environment: explicit setup from the
      DataCommunicator's rank/size (addr/port must then be reachable by all
      ranks).

    Args:
        data_communicator: Kratos DataCommunicator; defaults to the parallel
            environment's default one.
        addr, port: torch.distributed rendezvous address for the explicit
            setup paths.
        backend: torch.distributed backend; defaults to "gloo" (CPU-safe);
            pass "nccl" for GPU jobs.

    Returns:
        The initialized physicsnemo DistributedManager instance.
    """
    DistributedManager = _TryImportDistributedManager()

    if data_communicator is None:
        data_communicator = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()

    if DistributedManager.is_initialized():
        manager = DistributedManager()
        _CheckConsistency(manager, data_communicator)
        return manager

    if data_communicator.Size() == 1:
        DistributedManager.setup(rank=0, world_size=1, local_rank=0,
                                 addr=addr, port=port, backend=backend or "gloo")
        return DistributedManager()

    if any(name in os.environ for name in _LAUNCHER_ENV_VARIABLES):
        DistributedManager.initialize()
    else:
        import torch  # physicsnemo guarantees torch is present
        device_count = torch.cuda.device_count()
        local_rank = data_communicator.Rank() % device_count if device_count > 0 else 0
        DistributedManager.setup(
            rank=data_communicator.Rank(),
            world_size=data_communicator.Size(),
            local_rank=local_rank,
            addr=addr, port=port, backend=backend or "gloo")

    manager = DistributedManager()
    _CheckConsistency(manager, data_communicator)
    return manager


# --- MPI-aware field gathering (pure Kratos, no torch/physicsnemo) ----------


def _GetLocalContainer(model_part: Kratos.ModelPart, data_location: str):
    """Ghost-free local container for a data location (each entity is owned
    by exactly one rank, so gathering local containers never duplicates)."""
    local_mesh = model_part.GetCommunicator().LocalMesh()
    if data_location.startswith("node"):
        return local_mesh.Nodes
    if data_location.startswith("element"):
        return local_mesh.Elements
    if data_location.startswith("condition"):
        return local_mesh.Conditions
    raise RuntimeError(f"Unsupported data location \"{data_location}\".")


def GatherFieldToRank0(model_part: Kratos.ModelPart,
                       variable_name: str,
                       data_location: str,
                       data_communicator=None):
    """Gathers a field of a (possibly distributed) model part onto rank 0.

    Ghost entities are excluded (communicator LocalMesh containers), the
    per-rank blocks travel through DataCommunicator.Gatherv, and rank 0
    sorts the concatenated result by global entity id - the same order a
    serial run's containers have, so downstream consumers are
    rank-count-agnostic.

    Args:
        model_part: The model part (serial parts pass straight through).
        variable_name: Name of the Kratos variable.
        data_location: One of the tensor-adaptor data locations.
        data_communicator: Defaults to the model part communicator's one.

    Returns:
        (entity_ids, values) on rank 0 - values shaped
        (n_global_entities, *per_entity_shape) float64, ids int64 ascending;
        (None, None) on every other rank.
    """
    import numpy
    from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import (
        GetContainerTensorAdaptor)

    if data_communicator is None:
        data_communicator = model_part.GetCommunicator().GetDataCommunicator()

    container = _GetLocalContainer(model_part, data_location)
    variable = Kratos.KratosGlobals.GetVariable(variable_name)
    entity_ids = numpy.fromiter(
        (entity.Id for entity in container), dtype=numpy.int64, count=len(container))
    if len(container) > 0:
        values = numpy.array(
            GetContainerTensorAdaptor(
                container, data_location, variable, model_part.ProcessInfo).data,
            dtype=numpy.float64)
        trailing_shape = list(values.shape[1:])
    else:
        values = numpy.zeros((0,))
        trailing_shape = []

    if not data_communicator.IsDistributed():
        order = numpy.argsort(entity_ids, kind="stable")
        return entity_ids[order], values[order]

    gathered_shapes = data_communicator.GathervInts([int(n) for n in trailing_shape], 0)
    gathered_ids = data_communicator.GathervInts([int(i) for i in entity_ids], 0)
    gathered_values = data_communicator.GathervDoubles(
        [float(v) for v in values.ravel()], 0)

    if data_communicator.Rank() != 0:
        return None, None

    trailing = ()
    for rank, rank_ids in enumerate(gathered_ids):
        if len(rank_ids) > 0:
            trailing = tuple(gathered_shapes[rank])
            break

    all_ids = numpy.concatenate(
        [numpy.asarray(rank_ids, dtype=numpy.int64) for rank_ids in gathered_ids])
    all_values = numpy.concatenate(
        [numpy.asarray(rank_values, dtype=numpy.float64) for rank_values in gathered_values])
    all_values = all_values.reshape((len(all_ids),) + trailing)
    order = numpy.argsort(all_ids, kind="stable")
    return all_ids[order], all_values[order]


# --- MPI-aware mesh topology gathering (pure Kratos, no torch/physicsnemo) --

# Generic registered names used to reconstruct gathered entities on rank 0,
# keyed by geometry type. CAREFUL: "Element3D6N" is the PRISM; 3D surface
# geometries need the SurfaceElement3D*N names.
_ELEMENT_NAME_BY_GEOMETRY = {
    Kratos.GeometryData.KratosGeometryType.Kratos_Triangle2D3: "Element2D3N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Triangle2D6: "Element2D6N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Triangle3D3: "Element3D3N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Triangle3D6: "SurfaceElement3D6N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Quadrilateral2D4: "Element2D4N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Quadrilateral2D8: "Element2D8N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Quadrilateral2D9: "Element2D9N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Quadrilateral3D4: "SurfaceElement3D4N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Quadrilateral3D8: "SurfaceElement3D8N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Quadrilateral3D9: "SurfaceElement3D9N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Tetrahedra3D4: "Element3D4N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Tetrahedra3D10: "Element3D10N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Hexahedra3D8: "Element3D8N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Hexahedra3D20: "Element3D20N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Hexahedra3D27: "Element3D27N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Prism3D6: "Element3D6N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Prism3D15: "Element3D15N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Pyramid3D5: "Element3D5N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Pyramid3D13: "Element3D13N",
}
_CONDITION_NAME_BY_GEOMETRY = {
    Kratos.GeometryData.KratosGeometryType.Kratos_Triangle3D3: "SurfaceCondition3D3N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Triangle3D6: "SurfaceCondition3D6N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Quadrilateral3D4: "SurfaceCondition3D4N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Quadrilateral3D8: "SurfaceCondition3D8N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Quadrilateral3D9: "SurfaceCondition3D9N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Quadrilateral2D4: "PrismCondition2D4N",
    Kratos.GeometryData.KratosGeometryType.Kratos_Prism3D6: "PrismCondition3D6N",
}


def _GetEntityCreationName(geometry_type, source_container: str) -> str:
    table = _ELEMENT_NAME_BY_GEOMETRY if source_container == "Elements" else _CONDITION_NAME_BY_GEOMETRY
    name = table.get(geometry_type)
    if name is None:
        raise RuntimeError(
            f"No generic registered {source_container[:-1].lower()} name for geometry type "
            f"{geometry_type}; supported types: {sorted(str(t) for t in table)}.")
    return name


def _CheckSourceContainer(source_container):
    if source_container not in ("Elements", "Conditions", None):
        raise ValueError(
            f"Unsupported source container \"{source_container}\". "
            "Use \"Elements\", \"Conditions\" or None (nodes only).")


class GatheredMesh:
    """Topology of a (possibly distributed) model part, gathered on rank 0.

    Attributes (all id-sorted, ascending):
        node_ids: (N,) int64.        coordinates: (N, 3) float64.
        entity_ids: (E,) int64.      geometry_codes: (E,) int64
            (int(KratosGeometryType), round-trips through the enum).
        connectivity: list of E int64 arrays of node ids (full node lists,
            higher-order nodes included).
    """

    def __init__(self, node_ids, coordinates, entity_ids, geometry_codes, connectivity):
        self.node_ids = node_ids
        self.coordinates = coordinates
        self.entity_ids = entity_ids
        self.geometry_codes = geometry_codes
        self.connectivity = connectivity


def _SortMeshArrays(node_ids, coordinates, entity_ids, geometry_codes, counts, flat_connectivity):
    import numpy
    node_order = numpy.argsort(node_ids, kind="stable")
    offsets = numpy.zeros(len(counts) + 1, dtype=numpy.int64)
    numpy.cumsum(counts, out=offsets[1:])
    connectivity = [flat_connectivity[offsets[i]:offsets[i + 1]] for i in range(len(counts))]
    entity_order = numpy.argsort(entity_ids, kind="stable")
    return GatheredMesh(
        node_ids[node_order], coordinates[node_order],
        entity_ids[entity_order], geometry_codes[entity_order],
        [connectivity[int(i)] for i in entity_order])


def GatherMeshToRank0(model_part: Kratos.ModelPart,
                      source_container: str = "Elements",
                      data_communicator=None):
    """Gathers the ghost-free topology of a model part onto rank 0.

    Nodes come from the communicator LocalMesh (each node owned by exactly
    one rank); entities carry their geometry type and full node-id
    connectivity through Gatherv (counts + flat ids). Rank 0 returns a
    GatheredMesh with everything sorted by global id - the order a serial
    run's containers have. Other ranks return None. Serial parts pass
    through locally (same GatheredMesh contract).

    Args:
        model_part: The model part.
        source_container: "Elements" (default), "Conditions", or None for a
            nodes-only gather (empty entity arrays).
        data_communicator: Defaults to the model part communicator's one.
    """
    import numpy

    _CheckSourceContainer(source_container)
    if data_communicator is None:
        data_communicator = model_part.GetCommunicator().GetDataCommunicator()

    local_mesh = model_part.GetCommunicator().LocalMesh()
    nodes = local_mesh.Nodes
    node_ids = numpy.fromiter((node.Id for node in nodes), dtype=numpy.int64, count=len(nodes))
    if len(nodes) > 0:
        position_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(nodes, Kratos.Configuration.Current)
        position_ta.CollectData()
        coordinates = numpy.array(position_ta.data, dtype=numpy.float64).reshape(len(nodes), 3)
    else:
        coordinates = numpy.zeros((0, 3))

    entity_ids, geometry_codes, counts, flat_connectivity = [], [], [], []
    if source_container is not None:
        container = local_mesh.Elements if source_container == "Elements" else local_mesh.Conditions
        for entity in container:
            geometry = entity.GetGeometry()
            entity_ids.append(entity.Id)
            geometry_codes.append(int(geometry.GetGeometryType()))
            counts.append(len(geometry))
            flat_connectivity.extend(node.Id for node in geometry)
    entity_ids = numpy.asarray(entity_ids, dtype=numpy.int64)
    geometry_codes = numpy.asarray(geometry_codes, dtype=numpy.int64)
    counts = numpy.asarray(counts, dtype=numpy.int64)
    flat_connectivity = numpy.asarray(flat_connectivity, dtype=numpy.int64)

    if not data_communicator.IsDistributed():
        return _SortMeshArrays(node_ids, coordinates, entity_ids, geometry_codes, counts, flat_connectivity)

    gathered = [
        data_communicator.GathervInts([int(i) for i in node_ids], 0),
        data_communicator.GathervDoubles([float(v) for v in coordinates.ravel()], 0),
        data_communicator.GathervInts([int(i) for i in entity_ids], 0),
        data_communicator.GathervInts([int(i) for i in geometry_codes], 0),
        data_communicator.GathervInts([int(i) for i in counts], 0),
        data_communicator.GathervInts([int(i) for i in flat_connectivity], 0),
    ]
    if data_communicator.Rank() != 0:
        return None

    def concatenate(per_rank, dtype):
        return numpy.concatenate([numpy.asarray(block, dtype=dtype) for block in per_rank]) \
            if per_rank else numpy.zeros(0, dtype=dtype)

    all_node_ids = concatenate(gathered[0], numpy.int64)
    all_coordinates = concatenate(gathered[1], numpy.float64).reshape(len(all_node_ids), 3)
    return _SortMeshArrays(
        all_node_ids, all_coordinates,
        concatenate(gathered[2], numpy.int64), concatenate(gathered[3], numpy.int64),
        concatenate(gathered[4], numpy.int64), concatenate(gathered[5], numpy.int64))


class GatheredModelPart:
    """Result of GatherModelPartToRank0.

    Attributes:
        model: The Kratos.Model OWNING the rank-0 shadow part (keep it alive
            as long as the part is used); None off rank 0 and in serial runs.
        model_part: The shadow part on rank 0, the ORIGINAL part in serial
            runs, None on other ranks of a distributed run.
        field_specs: The effective (variable_name, data_location) pairs valid
            on model_part - Gauss-point locations are gathered, collapsed to
            their per-entity mean and re-labeled "element"/"condition"
            (matching what the serial mesh-bridge collection does).
    """

    def __init__(self, model, model_part, field_specs):
        self.model = model
        self.model_part = model_part
        self.field_specs = field_specs


def GatherModelPartToRank0(model_part: Kratos.ModelPart,
                           field_specs=(),
                           source_container: str = "Elements",
                           data_communicator=None) -> GatheredModelPart:
    """Reconstructs a distributed model part as a serial "shadow" on rank 0.

    The general primitive for running serial-only machinery (tessellation,
    exports, ...) on distributed data: topology through GatherMeshToRank0,
    fields through GatherFieldToRank0, entities recreated with generic
    registered names per geometry type. The shadow's containers are built
    id-ascending, reproducing a serial run's (id-sorted) container order,
    so anything computed from it is rank-count-agnostic.

    Collective on all ranks of the communicator. Serial parts pass through
    (the original part is returned, specs unchanged).

    Args:
        model_part: The model part to gather.
        field_specs: iterable of (variable_name, data_location) pairs to
            attach to the shadow.
        source_container: "Elements" (default) or "Conditions".
        data_communicator: Defaults to the model part communicator's one.
    """
    import numpy
    from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import (
        GetTensorAdaptor)

    if source_container not in ("Elements", "Conditions"):
        raise ValueError(
            f"Unsupported source container \"{source_container}\". Use \"Elements\" or \"Conditions\".")
    field_specs = [(str(name), str(location)) for name, location in field_specs]
    if data_communicator is None:
        data_communicator = model_part.GetCommunicator().GetDataCommunicator()

    if not data_communicator.IsDistributed():
        return GatheredModelPart(None, model_part, list(field_specs))

    mesh = GatherMeshToRank0(model_part, source_container, data_communicator)

    effective_specs = []
    gathered_values = []
    for variable_name, data_location in field_specs:
        _, values = GatherFieldToRank0(model_part, variable_name, data_location, data_communicator)
        if data_location.endswith("_gauss_point"):
            # Gauss points have no counterpart on the shadow: collapse to the
            # per-entity mean, exactly like the serial mesh-bridge collection.
            effective_location = data_location[:-len("_gauss_point")]
            if values is not None:
                values = values.mean(axis=1)
        else:
            effective_location = data_location
        effective_specs.append((variable_name, effective_location))
        gathered_values.append(values)

    if mesh is None:  # non-writing rank; all collectives are done
        return GatheredModelPart(None, None, effective_specs)

    shadow_model = Kratos.Model()
    shadow = shadow_model.CreateModelPart(model_part.Name)
    shadow.ProcessInfo[Kratos.DOMAIN_SIZE] = model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
    shadow.ProcessInfo[Kratos.STEP] = model_part.ProcessInfo[Kratos.STEP]
    shadow.ProcessInfo[Kratos.TIME] = model_part.ProcessInfo[Kratos.TIME]
    for variable_name, effective_location in effective_specs:
        if effective_location == "node_historical":
            shadow.AddNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(variable_name))
    shadow.SetBufferSize(1)

    for node_id, xyz in zip(mesh.node_ids, mesh.coordinates):
        shadow.CreateNewNode(int(node_id), float(xyz[0]), float(xyz[1]), float(xyz[2]))
    properties = shadow.CreateNewProperties(0)
    for entity_id, code, connectivity in zip(mesh.entity_ids, mesh.geometry_codes, mesh.connectivity):
        name = _GetEntityCreationName(
            Kratos.GeometryData.KratosGeometryType(int(code)), source_container)
        node_list = [int(i) for i in connectivity]
        if source_container == "Elements":
            shadow.CreateNewElement(name, int(entity_id), node_list, properties)
        else:
            shadow.CreateNewCondition(name, int(entity_id), node_list, properties)

    for (variable_name, effective_location), values in zip(effective_specs, gathered_values):
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        tensor_adaptor = GetTensorAdaptor(shadow, effective_location, variable)
        tensor_adaptor.data[:] = numpy.asarray(values).reshape(tensor_adaptor.data.shape).astype(
            tensor_adaptor.data.dtype, copy=False)
        tensor_adaptor.StoreData()

    return GatheredModelPart(shadow_model, shadow, effective_specs)


# --- Shared process groups and device meshes (torch/physicsnemo, lazy) ------


def _UnbindGlooDefaultGroupDevice() -> None:
    """physicsnemo's setup passes device_id (cuda on CUDA-visible machines)
    to torch's init_process_group, binding it to the default group; with the
    gloo backend that later breaks new_group (ProcessGroupGloo has no
    perform_nocolor_split). Unbinding restores the plain group-creation
    path; harmless when nothing is bound."""
    import torch.distributed as dist
    if not dist.is_initialized() or dist.get_backend() != "gloo":
        return
    try:
        group = dist.distributed_c10d._get_default_group()
        if getattr(group, "bound_device_id", None) is not None:
            group.bound_device_id = None
    except Exception:  # noqa: BLE001 - private torch surface; never fatal
        pass


def _GetInitializedManager():
    DistributedManager = _TryImportDistributedManager()
    if not DistributedManager.is_initialized():
        raise RuntimeError(
            "physicsnemo's DistributedManager is not initialized. Call "
            "distributed_utils.InitializeDistributedManager() first so torch and Kratos "
            "agree on ranks before creating process groups or device meshes.")
    _UnbindGlooDefaultGroupDevice()
    return DistributedManager()


def _CheckGroupConsistency(manager, group_name: str, sub_data_communicator) -> None:
    """Raises when the physicsnemo subgroup and the Kratos sub-communicator
    disagree about this rank's place in the group."""
    if (manager.group_rank(group_name) != sub_data_communicator.Rank()
            or manager.group_size(group_name) != sub_data_communicator.Size()):
        raise RuntimeError(
            f"Process group \"{group_name}\" is inconsistent between the two worlds: torch "
            f"sees rank {manager.group_rank(group_name)} of {manager.group_size(group_name)}, "
            f"Kratos sees rank {sub_data_communicator.Rank()} of {sub_data_communicator.Size()}. "
            "The group memberships diverged - check the parent communicators.")


def CreateMatchedProcessGroup(name: str,
                              size: int,
                              data_communicator=None,
                              parent_group_name=None,
                              verbose: bool = False):
    """Creates a physicsnemo process subgroup AND the matching registered
    Kratos sub-communicator over the same ranks.

    The torch side goes through DistributedManager.create_process_subgroup
    (collective; size must divide the parent group size; groups are blocks
    of consecutive ranks). The Kratos side splits the DataCommunicator with
    the same membership and registers it under the same name
    (ParallelEnvironment.GetDataCommunicator(name) retrieves it later). A
    consistency check makes any divergence loud.

    On a serial parent communicator only the torch-side subgroup is created
    and the parent communicator is returned unchanged.

    Args:
        name: Group name (must be new; also the Kratos registration name).
        size: Ranks per group (must divide the parent group size).
        data_communicator: Kratos parent; defaults to the parallel
            environment's default communicator.
        parent_group_name: Optional physicsnemo parent group to subdivide.
        verbose: Forwarded to physicsnemo.

    Returns:
        The Kratos sub-communicator this rank belongs to (the parent
        communicator itself on serial runs).
    """
    manager = _GetInitializedManager()
    DistributedManager = type(manager)
    if data_communicator is None:
        data_communicator = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()

    if parent_group_name is None:
        _CheckConsistency(manager, data_communicator)
    DistributedManager.create_process_subgroup(name, int(size), group_name=parent_group_name,
                                               verbose=verbose)

    if not data_communicator.IsDistributed():
        return data_communicator

    # Split the Kratos communicator with the subgroup's membership. The rank
    # lists live in physicsnemo's (private) _group_ranks; fall back to the
    # documented consecutive-block layout if that attribute ever disappears.
    from KratosMultiphysics.mpi import DataCommunicatorFactory
    rank = data_communicator.Rank()
    group_ranks = getattr(manager, "_group_ranks", {}).get(name)
    if group_ranks:
        color = next(i for i, ranks in enumerate(group_ranks) if rank in ranks)
        key = group_ranks[color].index(rank)
    else:  # consecutive blocks of `size` ranks
        color = rank // int(size)
        key = rank % int(size)
    sub_data_communicator = DataCommunicatorFactory.SplitAndRegister(
        data_communicator, int(color), int(key), name)

    _CheckGroupConsistency(manager, name, sub_data_communicator)
    return sub_data_communicator


def CreateMatchedProcessGroups(settings: Kratos.Parameters, data_communicator=None) -> dict:
    """Parameters-driven creation of matched process groups.

    Settings:
        {
            "process_groups" : [
                { "name": "model_parallel", "size": 2 },
                { "name": "mp_inner", "size": 1, "parent": "model_parallel" }
            ]
        }

    Entries are processed in order; "parent" names a previously created
    physicsnemo group to subdivide (empty = the world group).

    Returns:
        {name: Kratos sub-communicator} for every created group.
    """
    default_settings = Kratos.Parameters("""{
        "process_groups" : []
    }""")
    settings.ValidateAndAssignDefaults(default_settings)
    entry_defaults = Kratos.Parameters("""{
        "name"   : "PLEASE_SPECIFY_GROUP_NAME",
        "size"   : 1,
        "parent" : ""
    }""")

    groups = {}
    for i in range(settings["process_groups"].size()):
        entry = settings["process_groups"][i]
        entry.ValidateAndAssignDefaults(entry_defaults)
        parent = entry["parent"].GetString() or None
        groups[entry["name"].GetString()] = CreateMatchedProcessGroup(
            entry["name"].GetString(), entry["size"].GetInt(),
            data_communicator=data_communicator, parent_group_name=parent)
    return groups


def InitializeDeviceMesh(mesh_shape, mesh_dim_names, data_communicator=None):
    """Initializes physicsnemo's global device mesh, validated against Kratos.

    Args:
        mesh_shape: Ranks per mesh dimension; one entry may be -1 (inferred).
            The product must equal the communicator size.
        mesh_dim_names: One name per dimension (e.g. ("data", "model")).
        data_communicator: Defaults to the parallel environment's default.

    Returns:
        The torch.distributed.DeviceMesh.
    """
    manager = _GetInitializedManager()
    if data_communicator is None:
        data_communicator = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()

    mesh_shape = tuple(int(n) for n in mesh_shape)
    mesh_dim_names = tuple(str(n) for n in mesh_dim_names)
    if len(mesh_shape) != len(mesh_dim_names):
        raise ValueError(
            f"mesh_shape has {len(mesh_shape)} entries but mesh_dim_names has "
            f"{len(mesh_dim_names)}.")
    if sum(1 for n in mesh_shape if n == -1) > 1:
        raise ValueError(f"At most one mesh_shape entry may be -1, got {mesh_shape}.")

    world_size = data_communicator.Size()
    fixed_product = 1
    for n in mesh_shape:
        if n != -1:
            if n < 1:
                raise ValueError(f"mesh_shape entries must be positive (or one -1), got {mesh_shape}.")
            fixed_product *= n
    if -1 in mesh_shape:
        if world_size % fixed_product != 0:
            raise ValueError(
                f"mesh_shape {mesh_shape} cannot tile the {world_size} Kratos rank(s): "
                f"{world_size} is not divisible by {fixed_product}.")
    elif fixed_product != world_size:
        raise ValueError(
            f"mesh_shape {mesh_shape} implies {fixed_product} rank(s) but the Kratos "
            f"communicator has {world_size}.")

    _CheckConsistency(manager, data_communicator)
    return manager.initialize_mesh(mesh_shape, mesh_dim_names)
