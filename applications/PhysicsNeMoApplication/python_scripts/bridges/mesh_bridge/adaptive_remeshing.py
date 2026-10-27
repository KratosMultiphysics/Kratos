"""Surrogate-error-driven adaptive remeshing.

Closes the loop between error scoring and mesh adaptation:

- Volume path (Kratos/MMG): a per-node error - typically the assembled PDE
  residual of an ML-predicted field from solver_residuals.ResidualEvaluator
  - becomes a per-node target edge length (ComputeTargetSizeField), which
  drives MeshingApplication's MMG remesher through its scalar metric
  (RunMmgAdaptation). MeshingApplication with MMG support is an optional,
  lazily-imported dependency.
- Surface path (physicsnemo, pure torch): partition_cells clusters a
  triangle surface mesh around seed points; WeightedSurfacePartition samples
  the seeds with probability proportional to a per-cell error weight, so
  cluster density follows the error - the building block for error-adapted
  surface remeshing. RemeshSurface wraps physicsnemo.mesh.remeshing.remesh
  (isotropic clustering, Warp-backed as of 2.2, triangle surfaces only;
  its count targets output VERTICES and it drops cell data by design).

Kratos MMG facts this module encodes: the remesher reads the per-node
metric from METRIC_SCALAR (target edge length) unless the FIRST node
carries METRIC_TENSOR_<dim>D - so this module never calls MetricFastInit
(which would seed zero tensors everywhere and silently switch the mode) -
and it requires NODAL_H, computed here with FindNodalHNonHistoricalProcess.
"""

import numpy

import KratosMultiphysics as Kratos


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.adaptive_remeshing's surface utilities require torch, "
            "which could not be imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportPartitionCells():
    try:
        from physicsnemo.mesh.remeshing import partition_cells
        return partition_cells
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.adaptive_remeshing's surface utilities require "
            "physicsnemo, which could not be imported. Install it with e.g. "
            "'pip install nvidia-physicsnemo'.") from e


def _TryImportMeshingApplication():
    try:
        import KratosMultiphysics.MeshingApplication as MeshingApplication
    except ImportError as e:
        raise ImportError(
            "RunMmgAdaptation requires MeshingApplication, which is not available in "
            "this Kratos build. Compile it (add_app .../MeshingApplication) with MMG "
            "support (INCLUDE_MMG).") from e
    if not hasattr(MeshingApplication, "MmgProcess2D"):
        raise ImportError(
            "RunMmgAdaptation requires MeshingApplication compiled with MMG support "
            "(INCLUDE_MMG=ON); the available build exposes no MmgProcess.")
    return MeshingApplication


def NodalErrorArray(model_part: Kratos.ModelPart, nodal_residuals: dict) -> numpy.ndarray:
    """Collapses ResidualEvaluator.ComputeNodalResiduals output to one value
    per node (max over that node's DOFs), in model_part.Nodes iteration order.
    """
    per_node = {}
    for (node_id, _), value in nodal_residuals.items():
        per_node[node_id] = max(per_node.get(node_id, 0.0), value)
    return numpy.array([per_node.get(node.Id, 0.0) for node in model_part.Nodes])


def ComputeTargetSizeField(model_part: Kratos.ModelPart, nodal_error,
                           settings: Kratos.Parameters) -> numpy.ndarray:
    """Per-node target edge length from a per-node error indicator.

    The classic equidistribution rule: starting from the current local size
    NODAL_H (computed here), each node's target size is

        h_target = clip(h_current * (target_error / max(error, eps))**exponent,
                        minimal_size, maximal_size)

    so above-target errors refine (h shrinks) and below-target errors
    coarsen, with the exponent controlling how aggressively.

    Args:
        model_part: The model part to size (NODAL_H is (re)computed on it).
        nodal_error: (N,) array in model_part.Nodes iteration order, or the
            dict ComputeNodalResiduals returns (collapsed via NodalErrorArray).
        settings: Kratos Parameters:
            { "target_error": 1e-3, "exponent": 0.5,
              "minimal_size": 1e-4, "maximal_size": 1.0 }

    Returns:
        (N,) float64 array of target edge lengths.
    """
    defaults = Kratos.Parameters("""{
        "target_error" : 1e-3,
        "exponent"     : 0.5,
        "minimal_size" : 1e-4,
        "maximal_size" : 1.0
    }""")
    settings.ValidateAndAssignDefaults(defaults)
    target_error = settings["target_error"].GetDouble()
    exponent = settings["exponent"].GetDouble()
    minimal_size = settings["minimal_size"].GetDouble()
    maximal_size = settings["maximal_size"].GetDouble()
    KRATOS_EPS = 1e-30

    if isinstance(nodal_error, dict):
        nodal_error = NodalErrorArray(model_part, nodal_error)
    nodal_error = numpy.asarray(nodal_error, dtype=float)
    if nodal_error.shape != (model_part.NumberOfNodes(),):
        raise ValueError(
            f"nodal_error has shape {nodal_error.shape} but the model part has "
            f"{model_part.NumberOfNodes()} nodes.")

    Kratos.FindNodalHNonHistoricalProcess(model_part).Execute()
    current_h = numpy.array([node.GetValue(Kratos.NODAL_H) for node in model_part.Nodes])

    ratio = target_error / numpy.maximum(nodal_error, KRATOS_EPS)
    return numpy.clip(current_h * ratio**exponent, minimal_size, maximal_size)


def RunMmgAdaptation(model_part: Kratos.ModelPart, size_field,
                     mmg_parameters: Kratos.Parameters = None) -> None:
    """Remeshes a model part with MMG driven by a per-node scalar size field.

    Writes the sizes to METRIC_SCALAR (target edge length; the tensor metric
    is deliberately left unset so MMG stays in scalar mode), refreshes
    NODAL_H, and executes MmgProcess2D/3D according to the model part's
    DOMAIN_SIZE. Nodal values (historical and non-historical) are
    interpolated onto the new mesh by MMG itself.

    Args:
        model_part: The model part to remesh in place (root model part of
            the mesh; DOMAIN_SIZE must be set in its ProcessInfo).
        size_field: (N,) target edge lengths in Nodes iteration order.
        mmg_parameters: Optional extra/override Kratos Parameters for
            MmgProcess (merged over this module's defaults).
    """
    MeshingApplication = _TryImportMeshingApplication()

    size_field = numpy.asarray(size_field, dtype=float)
    if size_field.shape != (model_part.NumberOfNodes(),):
        raise ValueError(
            f"size_field has shape {size_field.shape} but the model part has "
            f"{model_part.NumberOfNodes()} nodes.")

    domain_size = model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
    if domain_size not in (2, 3):
        raise ValueError(f"DOMAIN_SIZE must be 2 or 3, got {domain_size}.")

    parameters = Kratos.Parameters("""{
        "discretization_type"        : "Standard",
        "framework"                  : "Eulerian",
        "interpolate_nodal_values"   : true,
        "interpolate_non_historical" : true,
        "save_external_files"        : false,
        "echo_level"                 : 0
    }""")
    if mmg_parameters is not None:
        for key in mmg_parameters.keys():
            if parameters.Has(key):
                parameters.RemoveValue(key)
            parameters.AddValue(key, mmg_parameters[key])

    Kratos.FindNodalHNonHistoricalProcess(model_part).Execute()
    for node, size in zip(model_part.Nodes, size_field):
        node.SetValue(MeshingApplication.METRIC_SCALAR, float(size))

    process_type = (MeshingApplication.MmgProcess2D if domain_size == 2
                    else MeshingApplication.MmgProcess3D)
    process_type(model_part, parameters).Execute()

    Kratos.FindNodalHNonHistoricalProcess(model_part).Execute()
    model_part.Set(Kratos.MODIFIED, True)


def WeightedSurfacePartition(mesh, n_clusters: int, weights=None, seed: int = -1):
    """Partitions a surface mesh's cells around error-weighted seed points.

    Seeds are cell centroids sampled without replacement with probability
    proportional to `weights` (uniform when None), then handed to
    physicsnemo's partition_cells - so regions with higher weight receive
    proportionally more clusters. Pure torch; always available.

    Args:
        mesh: physicsnemo.mesh.Mesh (any manifold dimension).
        n_clusters: Number of clusters/seeds (must not exceed the cell count).
        weights: Optional (n_cells,) nonnegative per-cell weights.
        seed: RNG seed for the sampling; -1 leaves the global RNG alone.

    Returns:
        (partition, seeds): the physicsnemo CellPartition (assignments,
        cluster_areas, cluster_normals, cluster_centroids) and the (n, D)
        seed coordinates used.
    """
    torch = _TryImportTorch()
    partition_cells = _TryImportPartitionCells()

    n_cells = int(mesh.cells.shape[0])
    if not 0 < n_clusters <= n_cells:
        raise ValueError(f"n_clusters must be in [1, {n_cells}], got {n_clusters}.")

    centroids = mesh.points[mesh.cells].mean(dim=1)  # (n_cells, D)
    if weights is None:
        probabilities = torch.ones(n_cells, dtype=centroids.dtype, device=centroids.device)
    else:
        probabilities = torch.as_tensor(weights, dtype=centroids.dtype,
                                        device=centroids.device).reshape(-1)
        if probabilities.shape[0] != n_cells or (probabilities < 0).any():
            raise ValueError(
                f"weights must be {n_cells} nonnegative values, got shape "
                f"{tuple(probabilities.shape)}.")

    generator = None
    if seed >= 0:
        generator = torch.Generator(device=centroids.device).manual_seed(seed)
    chosen = torch.multinomial(probabilities, n_clusters, replacement=False,
                               generator=generator)
    seeds = centroids[chosen].contiguous()
    return partition_cells(mesh, seeds), seeds


def RemeshSurface(mesh, n_vertices: int, max_iterations: int = 4,
                  transfer_point_data=False, resolution_field=None):
    """Isotropic surface remesh via physicsnemo.mesh.remeshing.remesh.

    Triangle surfaces embedded in 3D only (upstream raises NotImplementedError
    otherwise). Backed by Warp since physicsnemo 2.2 - the earlier pyacvd
    dependency is gone, and `n_vertices` now targets the output VERTEX count
    (it counted cells before 2.2, so a config carried over from then produces
    a mesh roughly twice as fine as intended).

    Args:
        mesh: physicsnemo.mesh.Mesh triangle surface.
        n_vertices: Target number of output vertices.
        max_iterations: Clustering iterations (upstream default 4).
        transfer_point_data: False (default) drops all point data; pass True
            or a key selection to carry fields onto the new topology. Cell
            data is always dropped - the topology is new.
        resolution_field: Optional per-point target-size field for graded
            (non-uniform) output.

    Returns:
        The remeshed physicsnemo.mesh.Mesh.
    """
    try:
        from physicsnemo.mesh.remeshing import remesh
    except ImportError as e:
        raise ImportError(
            "RemeshSurface requires physicsnemo (>= 2.2, Warp-backed), which could not "
            "be imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e
    return remesh(mesh, n_vertices, max_iterations=max_iterations,
                  transfer_point_data=transfer_point_data,
                  resolution_field=resolution_field)
