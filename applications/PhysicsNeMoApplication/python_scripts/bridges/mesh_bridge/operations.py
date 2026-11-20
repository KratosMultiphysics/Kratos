"""physicsnemo.mesh operations on Kratos meshes: curvature, repair,
subdivision, smoothing, extrusion and moments.

These are the mesh-processing primitives the bridge had no use for until
now, each answering a question that comes up once a surrogate is trained on
real geometry:

- WriteCurvatureFields puts mean and Gaussian curvature on the nodes, as a
  geometric FEATURE next to the signed distance field: a model learning a
  boundary-layer quantity wants to know where the surface bends.
- RepairMesh runs upstream's duplicate/degenerate/orientation pass before a
  signed-distance query or an export, where a malformed surface gives wrong
  answers rather than errors.
- SubdivideMesh and SmoothMesh refine and relax a surface.
- ExtrudeModelPart sweeps a 2-D Kratos case into a 3-D one and materializes
  it as real Kratos entities, so a planar solve can feed a model that only
  speaks volumes.
- IntegrateMoment is the P0 quadrature moment of two cell fields, which
  gives areas, centroids and inertia-style integrals in one call.

Every function takes or returns physicsnemo meshes as the rest of the mesh
bridge does; the model-part-facing wrappers say so in their names.

torch and physicsnemo are imported lazily.
"""

import numpy

import KratosMultiphysics as Kratos

_SUBDIVISION_SCHEMES = ("loop", "butterfly", "linear")


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.operations requires torch, which could "
            "not be imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportMeshOperations():
    try:
        from physicsnemo.mesh.curvature import (
            gaussian_curvature_vertices, mean_curvature_vertices)
        from physicsnemo.mesh.repair import repair_mesh
        from physicsnemo.mesh.smoothing import smooth_laplacian
        from physicsnemo.mesh.subdivision import (
            subdivide_butterfly, subdivide_linear, subdivide_loop)
        from physicsnemo.mesh.projections import extrude
        from physicsnemo.mesh.calculus.integration import integrate_moment
    except ImportError as e:
        raise ImportError(
            "The mesh operations require physicsnemo >= 2.2, which could not be "
            "imported. Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e
    return {
        "mean_curvature": mean_curvature_vertices,
        "gaussian_curvature": gaussian_curvature_vertices,
        "repair": repair_mesh,
        "smooth": smooth_laplacian,
        "loop": subdivide_loop,
        "butterfly": subdivide_butterfly,
        "linear": subdivide_linear,
        "extrude": extrude,
        "integrate_moment": integrate_moment,
    }


def ComputeCurvature(mesh, kind: str = "mean", include_boundary: bool = False):
    """Per-vertex curvature of a triangular surface.

    Args:
        mesh: A physicsnemo surface Mesh (spatial.BoundarySurface output).
        kind: "mean" or "gaussian".
        include_boundary: Mean curvature only - whether open-boundary
            vertices are computed rather than left at zero.

    Returns:
        An (n_points,) torch tensor.
    """
    operations = _TryImportMeshOperations()
    if kind == "mean":
        return operations["mean_curvature"](mesh, include_boundary=include_boundary)
    if kind == "gaussian":
        return operations["gaussian_curvature"](mesh)
    raise ValueError(
        f"Unsupported curvature \"{kind}\". Use \"mean\" or \"gaussian\".")


def WriteCurvatureFields(model_part: Kratos.ModelPart, settings: Kratos.Parameters) -> None:
    """Writes surface curvature onto the model part's nodes.

    The boundary surface keeps EVERY point of the tessellated mesh, so a
    vertex index is a node row and the values land on the right nodes with
    no provenance lookup.

    TRAP: a node on NO boundary facet - every interior node of a volume
    mesh - has zero vertex area, and upstream divides by it, so its
    curvature comes back NaN rather than zero or an error. Writing that
    onto a Kratos variable would poison every gather downstream, so
    non-finite values are replaced by zero here: an interior node has no
    surface curvature, which is what zero says.

    Args:
        model_part: The model part.
        settings: Kratos Parameters; defaults:
            mean_curvature_variable ("" = do not write it),
            gaussian_curvature_variable (""),
            output_location ("node_non_historical"),
            source_container ("Elements"),
            include_boundary (false),
            tessellation_mode ("smallest_id_diagonal").
    """
    from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
        domain_mesh_builder, spatial)

    default_settings = Kratos.Parameters("""{
        "mean_curvature_variable"      : "",
        "gaussian_curvature_variable"  : "",
        "output_location"              : "node_non_historical",
        "source_container"             : "Elements",
        "include_boundary"             : false,
        "tessellation_mode"            : "smallest_id_diagonal"
    }""")
    settings.ValidateAndAssignDefaults(default_settings)

    mean_name = settings["mean_curvature_variable"].GetString()
    gaussian_name = settings["gaussian_curvature_variable"].GetString()
    if not mean_name and not gaussian_name:
        raise ValueError(
            "Neither \"mean_curvature_variable\" nor \"gaussian_curvature_variable\" "
            "is set, so there is nothing to write.")

    mesh, _ = domain_mesh_builder.BuildMesh(
        model_part, source_container=settings["source_container"].GetString(),
        tessellation_mode=settings["tessellation_mode"].GetString())
    surface = spatial.BoundarySurface(mesh)
    include_boundary = settings["include_boundary"].GetBool()
    location = settings["output_location"].GetString()

    for name, kind in ((mean_name, "mean"), (gaussian_name, "gaussian")):
        if not name:
            continue
        values = numpy.asarray(
            ComputeCurvature(surface, kind, include_boundary).detach().cpu().numpy(),
            dtype=numpy.float64)
        values = numpy.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
        _WriteNodalValues(model_part, name, location, values)


def _WriteNodalValues(model_part, variable_name, data_location, values) -> None:
    """Writes an (n_nodes,) array onto the nodes, historical or not."""
    variable = Kratos.KratosGlobals.GetVariable(variable_name)
    if len(values) != model_part.NumberOfNodes():
        raise ValueError(
            f"Got {len(values)} values for {model_part.NumberOfNodes()} nodes.")
    if data_location == "node_historical":
        for node, value in zip(model_part.Nodes, values):
            node.SetSolutionStepValue(variable, float(value))
    elif data_location == "node_non_historical":
        for node, value in zip(model_part.Nodes, values):
            node.SetValue(variable, float(value))
    else:
        raise ValueError(
            f"Unsupported data location \"{data_location}\"; curvature is nodal "
            "(\"node_historical\" or \"node_non_historical\").")


def RepairMesh(mesh, **options):
    """Upstream's mesh repair pass, with its report.

    Worth running before a signed-distance query or an export: duplicate
    points, degenerate cells and inconsistent orientation make those give
    WRONG ANSWERS rather than errors.

    Returns:
        (repaired_mesh, report) - the report counting what was changed.
    """
    operations = _TryImportMeshOperations()
    return operations["repair"](mesh, **options)


def SubdivideMesh(mesh, scheme: str = "loop", levels: int = 1):
    """Refines a triangular surface, quadrupling its cells per level.

    Schemes: "loop" and "butterfly" (smoothing subdivisions that move the
    points) and "linear" (pure refinement, points stay on the surface).
    """
    operations = _TryImportMeshOperations()
    if scheme not in _SUBDIVISION_SCHEMES:
        raise ValueError(
            f"Unsupported subdivision scheme \"{scheme}\". Supported: "
            f"{', '.join(_SUBDIVISION_SCHEMES)}.")
    if levels < 1:
        raise ValueError(f"\"levels\" must be >= 1 [ levels = {levels} ].")
    for _ in range(levels):
        mesh = operations[scheme](mesh)
    return mesh


def SmoothMesh(mesh, n_iter: int = 20, relaxation_factor: float = 0.01, **options):
    """Laplacian smoothing of a surface, boundaries preserved by default."""
    operations = _TryImportMeshOperations()
    return operations["smooth"](
        mesh, n_iter=n_iter, relaxation_factor=relaxation_factor, **options)


def IntegrateMoment(mesh, left, right, aligned_dims: int = 0):
    """The P0 quadrature moment sum_c |cell_c| left_c (x) right_c.

    With left = right = ones this is the mesh's measure (its area for a
    surface, its volume for a volume mesh), which is the cheapest way to
    check a generated or repaired mesh is what it claims to be.
    """
    operations = _TryImportMeshOperations()
    return operations["integrate_moment"](mesh, left, right, aligned_dims=aligned_dims)


def ExtrudeMesh(mesh, vector, capping: bool = False):
    """Sweeps a mesh along a vector, one dimension up."""
    operations = _TryImportMeshOperations()
    return operations["extrude"](
        mesh, vector=vector, capping=capping, allow_new_spatial_dims=True)


def ExtrudeModelPart(model: Kratos.Model, source_model_part: Kratos.ModelPart,
                     vector, new_model_part_name: str, capping: bool = False,
                     source_container: str = "Elements",
                     settings: Kratos.Parameters = None) -> Kratos.ModelPart:
    """Sweeps a 2-D Kratos case into a 3-D one, as real Kratos entities.

    The planar model part is tessellated, extruded along the vector, and
    materialized through generate.PopulateModelPartFromMesh - so the result
    is a model part a solver or another bridge can consume, not just an
    array. Upstream decomposes the sweep into SIMPLICES rather than prisms,
    so a swept triangle arrives as three tetrahedra, and the element name
    follows the produced cell width unless the settings override it.

    An in-plane vector produces a flat, zero-measure sweep; the measure is
    checked and refused rather than materialized as a degenerate part.

    Args:
        model: The Kratos Model the new part is created in.
        source_model_part: The planar part to sweep.
        vector: The extrusion vector, e.g. [0.0, 0.0, 1.0].
        new_model_part_name: Name for the created part; it must not exist.
        capping: Whether to close the swept volume's ends.
        source_container: Which container the source is tessellated from.
        settings: Forwarded to PopulateModelPartFromMesh (element_name,
            historical_variables, buffer_size, ...).

    Returns:
        The created Kratos.ModelPart.
    """
    from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
        domain_mesh_builder, generate)

    mesh, _ = domain_mesh_builder.BuildMesh(
        source_model_part, source_container=source_container)
    swept = ExtrudeMesh(mesh, vector, capping=capping)
    torch = _TryImportTorch()
    ones = torch.ones(int(swept.cells.shape[0]), dtype=torch.as_tensor(swept.points).dtype)
    measure = float(IntegrateMoment(swept, ones, ones))
    if abs(measure) < 1e-12:
        raise ValueError(
            f"The extrusion produced a mesh of measure {measure:g}: the vector "
            f"{list(vector)} lies in the source's own plane, so nothing was swept. "
            "Give a vector with a component out of that plane.")
    return generate.PopulateModelPartFromMesh(
        model, new_model_part_name, swept, settings)
