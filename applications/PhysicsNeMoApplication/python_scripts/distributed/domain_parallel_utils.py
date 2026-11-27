"""Domain parallelism over physicsnemo's ``ShardTensor``: ONE tensor split
across the Kratos ranks, each holding the rows (or the grid slab) it owns.

This is the other kind of parallelism. ``graph_partition_utils`` is DATA
parallel - every rank runs the whole model on its own subgraph and the
gradients are averaged. Here every rank holds a piece of one tensor, and
physicsnemo's registered handlers make the model's own ops mesh-aware: a
pointwise layer runs on the local rows, a convolution exchanges the halo it
needs from its neighbours, and a reduction (``mean``, ``sum``) becomes a
mesh-wide quantity whose backward already delivers the full gradient on
every rank. It is what one wants when a single Kratos mesh or grid does not
fit on one device.

What runs here, over CPU/gloo, with no ``DistributedManager`` at all:
construction from uneven per-rank rows (Kratos partitions are never
``torch.chunk``-shaped), ``full_tensor``/``redistribute``, autograd through
all of it, and the halo-based ``conv``/pool/interpolate/normalization/
linear/reduction/view handlers. What stays CUDA-only upstream, honestly:
``scatter_tensor`` (needs the DistributedManager), ``sharding_shapes="infer"``
(hard-codes ``device="cuda"`` in its shape exchange - the reason
``ShardLocalRows`` exchanges shapes itself), ring attention, the sharded
kNN/radius search and the NCCL transport. The MPI suite asserts the shipped
paths at np=2/3 over gloo, exactly the status the FSDP2 checkpoints have.

Two traps this module encodes:

- ``DistributedManager.initialize_mesh`` (hence ``distributed_utils.
  InitializeDeviceMesh``) builds a CUDA mesh whenever a GPU is visible, even
  under gloo - two ranks bound to one GPU. ``InitializeDomainMesh`` goes to
  ``torch.distributed.device_mesh.init_device_mesh`` directly.
- Upstream registers the ShardTensor op handlers only when
  ``torch.cuda.is_available()``; on a CUDA-less host they must be registered
  explicitly, which ``_TryImportShardTensor`` does.
- Do NOT all-reduce parameter gradients after a loss over a ShardTensor:
  the reduction handler's backward has already summed them (measured: every
  rank holds the serial gradient to 3e-8). A DDP-style average on top would
  be wrong by the rank count.

torch and physicsnemo are imported lazily.
"""

import numpy

import KratosMultiphysics as Kratos

from KratosMultiphysics.PhysicsNeMoApplication.distributed import graph_partition_utils
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.inference_process import (
    GatherInputFields, SynchronizeOutputFields, WriteOutputFields)

_OPS_REGISTERED = False


def _TryImportShardTensor():
    """(torch, ShardTensor), with the ShardTensor op handlers registered."""
    global _OPS_REGISTERED
    try:
        import torch
        import torch.distributed  # noqa: F401
        import physicsnemo.domain_parallel as domain_parallel
        from physicsnemo.domain_parallel import ShardTensor
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.domain_parallel_utils requires torch (>= 2.6) and "
            "physicsnemo with physicsnemo.domain_parallel, which could not be imported. "
            "Install them with e.g. 'pip install torch nvidia-physicsnemo'.") from e
    if not _OPS_REGISTERED:
        # upstream only auto-registers when a CUDA device is visible
        domain_parallel.register_custom_ops()
        _OPS_REGISTERED = True
    return torch, ShardTensor


def InitializeDomainMesh(data_communicator=None, device_type: str = "cpu",
                         mesh_dim_name: str = "domain", backend: str = "gloo",
                         port: str = "29500"):
    """A one-dimensional torch DeviceMesh over the Kratos ranks.

    Built with torch's ``init_device_mesh`` on top of
    ``graph_partition_utils.InitializeTorchProcessGroup`` (idempotent), NOT
    with physicsnemo's ``DistributedManager.initialize_mesh``: that one puts
    the mesh on CUDA whenever a GPU is visible, regardless of the backend -
    on a single-GPU box two gloo ranks end up bound to the same device.

    Build it ONCE per program: creating a mesh creates process sub-groups,
    and doing that repeatedly while other groups are live deadlocks.

    Args:
        data_communicator: Defaults to the parallel environment's default.
        device_type: "cpu" (the verified path) or "cuda".
        mesh_dim_name: The single mesh dimension's name.
        backend/port: Forwarded to InitializeTorchProcessGroup when no
            torch process group exists yet.

    Returns:
        The torch.distributed.DeviceMesh.
    """
    torch, _ = _TryImportShardTensor()
    from torch.distributed.device_mesh import init_device_mesh

    if device_type not in ("cpu", "cuda"):
        raise ValueError(f"device_type must be \"cpu\" or \"cuda\", got \"{device_type}\".")
    _, world_size = graph_partition_utils.InitializeTorchProcessGroup(
        data_communicator, backend=backend, port=port)
    return init_device_mesh(device_type, (world_size,), mesh_dim_names=(str(mesh_dim_name),))


def ShardLocalRows(local_rows, mesh, dim: int = 0):
    """A ShardTensor from each rank's OWN rows along ``dim``.

    Per-rank shapes may differ (a Kratos partition owns however many nodes
    it owns). They are exchanged with ``all_gather_object`` over the mesh's
    group and handed to ``from_local`` explicitly, because upstream's
    ``sharding_shapes="infer"`` does that exchange on CUDA unconditionally.

    Args:
        local_rows: This rank's (n_local, ...) numpy array or torch tensor.
        mesh: The DeviceMesh from InitializeDomainMesh.
        dim: The sharded axis.

    Returns:
        The ShardTensor with placements [Shard(dim)].
    """
    torch, ShardTensor = _TryImportShardTensor()
    import torch.distributed as distributed
    from torch.distributed.tensor import Shard

    local = torch.as_tensor(local_rows)
    if local.ndim == 0:
        raise ValueError("local_rows must have at least one axis to shard along.")
    dim = dim % local.ndim
    group = mesh.get_group()
    shapes = [None] * distributed.get_world_size(group)
    distributed.all_gather_object(shapes, tuple(int(n) for n in local.shape), group=group)
    for other in shapes:
        if len(other) != local.ndim or any(
                a != b for axis, (a, b) in enumerate(zip(other, local.shape)) if axis != dim):
            raise ValueError(
                f"Local shapes {shapes} cannot be concatenated along axis {dim}; every "
                "axis but the sharded one must agree across ranks.")
    # sharding_shapes is keyed by MESH dimension (0: the only one here),
    # not by the tensor axis - a {tensor_dim: ...} dict raises KeyError
    return ShardTensor.from_local(local, mesh, [Shard(dim)], sharding_shapes={0: shapes})


def ShardKratosField(model_part: Kratos.ModelPart, field_specs, mesh):
    """The OWNED rows of a model part's fields as one ShardTensor.

    Rows are gathered ghost-free (``GatherInputFields(local_only=True)``), so
    the ranks' shards partition the global entity set exactly and the global
    row order of the ShardTensor is rank-major. The owned entity ids come
    back with it, which is what lets ``full_tensor()`` be reordered into the
    id-sorted layout a serial run has.

    Args:
        model_part: A (possibly distributed) model part.
        field_specs: [(variable_name, data_location)], as InferenceProcess.
        mesh: The DeviceMesh from InitializeDomainMesh.

    Returns:
        (shard_tensor, owned_ids): the (n_owned, total_width) shard and the
        int64 ids of this rank's rows, in row order.
    """
    from KratosMultiphysics.PhysicsNeMoApplication.distributed.distributed_utils import (
        _GetLocalContainer)

    torch, _ = _TryImportShardTensor()
    field_specs = list(field_specs)
    if not field_specs:
        raise ValueError("field_specs must name at least one field.")
    inputs, _ = GatherInputFields(model_part, field_specs, local_only=True)
    features = torch.cat(inputs, dim=-1)
    container = _GetLocalContainer(model_part, field_specs[0][1])
    owned_ids = numpy.fromiter((entity.Id for entity in container), dtype=numpy.int64,
                               count=len(container))
    return ShardLocalRows(features, mesh, dim=0), owned_ids


def WriteShardedField(shard_tensor, model_part: Kratos.ModelPart, field_specs,
                      normalization=None, scale_only: bool = False) -> None:
    """Writes a ShardTensor's local rows onto this rank's OWNED entities and
    refreshes the ghosts from their owners.

    The rows must be the ones ShardKratosField produced (same entities, same
    order). ``normalization`` is a model_registry.LoadOutputNormalization
    entry, applied as everywhere else on the write path.
    """
    local = shard_tensor.to_local().detach()
    WriteOutputFields(model_part, list(field_specs), local, int(local.shape[0]),
                      local_only=True, normalization=normalization, scale_only=scale_only)
    SynchronizeOutputFields(model_part, list(field_specs))


def ShardGridAlongAxis(grid, mesh, axis: int):
    """A replicated (C, *spatial) grid as a ShardTensor split along ``axis``.

    Every rank passes the SAME full grid (e.g. the output of
    grid_bridge.SampleFieldsOnGrid on rank 0, broadcast) and keeps its
    ``torch.chunk`` piece - the layout physicsnemo's halo-exchanging
    convolution and pooling handlers reassemble exactly.

    Args:
        grid: The full grid, numpy or torch, identical on every rank.
        mesh: The DeviceMesh from InitializeDomainMesh.
        axis: The spatial axis to split (never 0, the channel axis).

    Returns:
        The ShardTensor with placements [Shard(axis)].
    """
    torch, ShardTensor = _TryImportShardTensor()
    from torch.distributed.tensor import Shard

    grid = torch.as_tensor(grid)
    axis = axis % grid.ndim
    if axis == 0:
        raise ValueError("axis 0 is the channel axis; shard a spatial axis.")
    world_size, rank = mesh.size(), mesh.get_local_rank()
    if grid.shape[axis] < world_size:
        raise ValueError(
            f"Axis {axis} has {grid.shape[axis]} points, fewer than the {world_size} rank(s) "
            "it would be split over.")
    pieces = torch.chunk(grid, world_size, dim=axis)
    return ShardTensor.from_local(pieces[rank].contiguous(), mesh, [Shard(axis)],
                                  sharding_shapes="chunk", global_shape=tuple(grid.shape))


def MeshWideValue(sharded_scalar) -> float:
    """The mesh-wide value of a reduced ShardTensor (a loss), as a float.

    ``float(loss)`` on a ShardTensor reads this rank's LOCAL partial, not the
    global reduction; ``full_tensor()`` is the collective that returns the
    value every rank agrees on. Every rank must call it.
    """
    value = sharded_scalar
    if hasattr(value, "full_tensor"):
        value = value.full_tensor()
    return float(value.detach().reshape(-1)[0]) if hasattr(value, "detach") else float(value)
