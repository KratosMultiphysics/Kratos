"""Mesh generation and repair from implicit geometry (physicsnemo >= 2.2).

The mesh bridge has always run one way - a Kratos ModelPart becomes a
physicsnemo Mesh, and predictions are scattered back onto entities that
already exist. This module adds the missing direction: meshes are
*generated* (from an implicit function or a level-set field) and then
*materialized as real Kratos entities*, so a geometry defined by an SDF can
be solved on.

It composes with what already ships:

    ModelPart --spatial.SampleSignedDistanceOnGrid--> level-set field
              --SurfaceFromLevelSet--------------> outward-wound surface
    phi       --GenerateImplicitDomain-----------> volume mesh (tets)
              --PopulateModelPartFromMesh--------> a solvable ModelPart
              --adaptive_remeshing.RunMmgAdaptation--> quality cleanup
    old part  --mapping_bridge.MappingBridge------> fields on the new mesh

Three upstream behaviours are handled here rather than passed through:

- **`mesh_implicit_domain` requires grad mode.** Under `torch.no_grad()` its
  boundary Newton projection degenerates and the coverage guard trips (or,
  with the guard disabled, it silently returns a worse mesh). Since any
  deployment process may have wrapped the world in `no_grad`,
  GenerateImplicitDomain forces `torch.enable_grad()` itself.
- **`marching_cubes` output is always CPU float32 and detached**, whatever
  you feed it - stated here because nothing upstream warns.
- **`fill_interior` is 2D-only.** A 3D triangle surface raises a bare
  NotImplementedError upstream (exact 3D boundary recovery is "planned").
  FillSurfaceWithTetrahedra tries upstream first anyway - so the day it
  lands this module inherits it - and otherwise falls back to a Delaunay
  tetrahedralization carved by the winding-number sign, validated against
  the input's own volume and boundary area.

torch/physicsnemo are optional runtime dependencies - imported lazily.
"""

import numpy

import KratosMultiphysics as Kratos


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.generate requires torch, which could not "
            "be imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportGenerators():
    try:
        from physicsnemo.mesh.generate import (
            marching_cubes, mesh_implicit_domain, refit_mesh_to_implicit)
        return mesh_implicit_domain, marching_cubes, refit_mesh_to_implicit
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.generate requires physicsnemo >= 2.2 "
            "(its mesh generators landed in 2.2), which could not be imported. Install "
            "it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def _TryImportFillInterior():
    try:
        from physicsnemo.mesh.tessellation import fill_interior
        return fill_interior
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.generate requires physicsnemo >= 2.2 "
            "(fill_interior landed in 2.2), which could not be imported. Install it "
            "with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def _TryImportTetgen():
    try:
        import tetgen
        return tetgen
    except ImportError as e:
        raise ImportError(
            "\"method\": \"tetgen\" requires the tetgen package (Python bindings of "
            "TetGen - note TetGen itself is AGPL-licensed, which is why this backend "
            "is an explicit opt-in), which could not be imported. Install it with "
            "e.g. 'pip install tetgen'.") from e


def SdfPrimitives():
    """The implicit-geometry building blocks, as a namespace dict.

    Returns physicsnemo's `sdf_sphere`, `sdf_box`, `sdf_polygon_2d` and the
    `sdf_union`/`sdf_intersection`/`sdf_difference` combinators. They are
    plain closures over points, negative inside, differentiable, and
    dtype/device polymorphic - so they compose with anything here.
    """
    try:
        from physicsnemo.mesh.generate import (
            sdf_box, sdf_difference, sdf_intersection, sdf_polygon_2d, sdf_sphere,
            sdf_union)
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.generate requires physicsnemo >= 2.2 "
            "(the sdf primitives landed in 2.2), which could not be imported. Install "
            "it with e.g. 'pip install -U nvidia-physicsnemo'.") from e
    return {"sphere": sdf_sphere, "box": sdf_box, "polygon_2d": sdf_polygon_2d,
            "union": sdf_union, "intersection": sdf_intersection,
            "difference": sdf_difference}


def GenerateImplicitDomain(phi, bounds, h, settings: Kratos.Parameters = None,
                           feature_points=None):
    """Generates a volume mesh filling the region where `phi` is negative.

    Args:
        phi: callable(points (..., D)) -> (...) signed values, negative
            inside. Must be autograd-differentiable (the boundary
            projection differentiates it).
        bounds: (low, high), each a length-D sequence bounding the domain.
            Memory scales with prod((high - low) / h), so keep it tight.
        h: Target edge length of the background lattice.
        settings: Optional Kratos Parameters:
            {
                "reconnect"    : "flips",   // or "none" (frozen topology)
                "iters"        : 60,
                "dtype"        : "float64", // float32 hurts the conditioning
                "seed"         : 0,
                "full_output"  : false
            }
        feature_points: Optional (P, D) points to pin exactly - needed for
            sharp corners, which otherwise trip the coverage guard.

    Returns:
        The physicsnemo Mesh, or (mesh, diagnostics) with "full_output".

    Notes:
        Runs under an explicit `torch.enable_grad()`: upstream needs grad
        mode even though nothing is being trained, and degrades silently
        without it. CPU generation is bit-exact reproducible; CUDA is NOT
        (atomics) and was slower at these sizes, so this stays on the CPU.
    """
    torch = _TryImportTorch()
    mesh_implicit_domain, _, _ = _TryImportGenerators()

    defaults = Kratos.Parameters("""{
        "reconnect"   : "flips",
        "iters"       : 60,
        "dtype"       : "float64",
        "seed"        : 0,
        "full_output" : false
    }""")
    settings = Kratos.Parameters("{}") if settings is None else settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)

    reconnect = settings["reconnect"].GetString()
    if reconnect not in ("flips", "none"):
        raise ValueError(
            f"Unknown reconnect strategy \"{reconnect}\". Use \"flips\" or \"none\".")
    dtype_name = settings["dtype"].GetString()
    dtypes = {"float64": torch.float64, "float32": torch.float32}
    if dtype_name not in dtypes:
        raise ValueError(
            f"Unsupported dtype \"{dtype_name}\". Use one of {tuple(dtypes)}.")

    low, high = (numpy.asarray(b, dtype=float).reshape(-1) for b in bounds)
    if low.shape != high.shape or numpy.any(high <= low):
        raise ValueError(
            f"bounds must be (low, high) of equal length with high > low, got "
            f"{low.tolist()} and {high.tolist()}.")

    keywords = {"reconnect": reconnect, "iters": settings["iters"].GetInt(),
                "dtype": dtypes[dtype_name], "device": "cpu",
                "seed": settings["seed"].GetInt(),
                "full_output": settings["full_output"].GetBool()}
    if feature_points is not None:
        keywords["feature_points"] = torch.as_tensor(
            numpy.asarray(feature_points, dtype=float), dtype=dtypes[dtype_name])

    # enable_grad, not no_grad: upstream's boundary projection differentiates
    # phi, and silently produces a worse mesh when grad mode is off
    with torch.enable_grad():
        return mesh_implicit_domain(phi, (low.tolist(), high.tolist()), float(h), **keywords)


def SurfaceFromLevelSet(field, bounding_box=None, threshold: float = 0.0):
    """Extracts an outward-wound triangle surface from a 3D level-set field.

    The natural partner of spatial.SampleSignedDistanceOnGrid: sample a
    model part's SDF onto a lattice, then pull the zero level set back out
    as a surface.

    Unlike `Mesh.get_boundary_mesh()` - whose triangles are wound
    inconsistently, giving a closed-surface signed volume of zero - marching
    cubes produces a consistently outward-oriented surface, so the result
    can be used directly for signed-distance queries.

    Args:
        field: (nx, ny, nz) scalar field, or (1, nx, ny, nz) as
            SampleSignedDistanceOnGrid returns.
        bounding_box: Optional (low, high) giving the physical extent; without
            it the surface comes back in grid-index space.
        threshold: The level to extract (0.0 for an SDF's surface).

    Returns:
        A physicsnemo Mesh of triangles. Always CPU float32 and detached -
        upstream converts regardless of the input's dtype or device.
    """
    torch = _TryImportTorch()
    _, marching_cubes, _ = _TryImportGenerators()

    field = torch.as_tensor(numpy.asarray(field))
    if field.dim() == 4 and field.shape[0] == 1:
        field = field[0]
    if field.dim() != 3:
        raise ValueError(
            f"marching cubes needs a 3D scalar field (nx, ny, nz), got shape "
            f"{tuple(field.shape)}.")

    keywords = {}
    if bounding_box is not None:
        low, high = (numpy.asarray(b, dtype=float).reshape(-1) for b in bounding_box)
        keywords["coords"] = tuple(
            torch.linspace(float(low[axis]), float(high[axis]), int(field.shape[axis]))
            for axis in range(3))
    return marching_cubes(field, float(threshold), **keywords)


def FillBoundaryLoop(boundary, settings: Kratos.Parameters = None):
    """Fills closed 2D boundary loops with quality-guaranteed triangles.

    Args:
        boundary: physicsnemo Mesh of EDGES (cells of shape (E, 2)) forming
            one or more closed loops; nesting (holes) is resolved upstream.
        settings: Optional Kratos Parameters:
            {
                "max_cell_size"     : 0.0,   // triangle AREA, 0 = unconstrained
                "min_angle_degrees" : 30.0,  // upstream caps this at 33
                "smooth_iterations" : 0
            }

    Returns:
        A physicsnemo Mesh of triangles meeting the minimum-angle guarantee.
        The input loop's vertices are preserved as the leading rows.
    """
    torch = _TryImportTorch()
    fill_interior = _TryImportFillInterior()

    defaults = Kratos.Parameters("""{
        "max_cell_size"     : 0.0,
        "min_angle_degrees" : 30.0,
        "smooth_iterations" : 0
    }""")
    settings = Kratos.Parameters("{}") if settings is None else settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)

    if int(boundary.cells.shape[1]) != 2:
        raise ValueError(
            "FillBoundaryLoop needs a 2D boundary of EDGES (cells of shape (E, 2)); got "
            f"cells of shape {tuple(boundary.cells.shape)}. To fill a 3D triangle "
            "surface with tetrahedra use FillSurfaceWithTetrahedra, or "
            "GenerateImplicitDomain when the geometry is an implicit function.")

    max_cell_size = settings["max_cell_size"].GetDouble()
    return fill_interior(
        boundary,
        max_cell_size=(max_cell_size if max_cell_size > 0.0 else None),
        min_angle_degrees=settings["min_angle_degrees"].GetDouble(),
        smooth_iterations=settings["smooth_iterations"].GetInt())


def _SurfaceMetrics(points, triangles):
    """(enclosed volume, area) of a closed, consistently wound triangle surface."""
    a, b, c = points[triangles[:, 0]], points[triangles[:, 1]], points[triangles[:, 2]]
    volume = float(numpy.einsum("ij,ij->i", a, numpy.cross(b, c)).sum() / 6.0)
    area = float(numpy.linalg.norm(numpy.cross(b - a, c - a), axis=1).sum() / 2.0)
    return volume, area


def _BoundaryFaces(tetrahedra):
    """The faces used by exactly one tetrahedron, as sorted index triples."""
    counts = {}
    for tet in tetrahedra:
        for corners in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            key = tuple(sorted(int(tet[i]) for i in corners))
            counts[key] = counts.get(key, 0) + 1
    return numpy.array([face for face, n in counts.items() if n == 1], dtype=numpy.int64)


def _IsEdgeManifold(faces):
    """True when every undirected edge of `faces` is shared by exactly two of them."""
    counts = {}
    for face in faces:
        for i, j in ((0, 1), (1, 2), (2, 0)):
            key = (int(face[i]), int(face[j])) if face[i] < face[j] else (int(face[j]), int(face[i]))
            counts[key] = counts.get(key, 0) + 1
    return bool(counts) and set(counts.values()) == {2}


def _CarveDelaunay(surface, points):
    """Delaunay tetrahedra whose centroids fall inside the surface.

    scipy tetrahedralizes the convex hull of the points; the winding-number
    sign of each centroid then discards whatever lies outside the actual
    (possibly non-convex) solid. The points are passed through untouched, so
    every input vertex survives bit-identically.
    """
    from scipy.spatial import Delaunay

    torch = _TryImportTorch()
    try:
        from physicsnemo.mesh.spatial import signed_distance_field
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.generate requires physicsnemo >= 2.2 "
            "(mesh.spatial.signed_distance_field landed in 2.2), which could not be "
            "imported. Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e

    candidates = numpy.asarray(Delaunay(points).simplices, dtype=numpy.int64)
    centroids = numpy.ascontiguousarray(points[candidates].mean(axis=1))
    # warp's from_torch rejects non-contiguous input, and the winding-number
    # sign is what makes this robust on non-convex solids.
    distances, _, _ = signed_distance_field(
        surface, torch.as_tensor(centroids).contiguous(), use_sign_winding_number=True)
    tetrahedra = candidates[numpy.asarray(distances.detach().cpu()) < 0.0]

    # Delaunay does not wind its simplices consistently; Kratos needs positive
    # Jacobians, so flip the negative ones. Flat slivers (coplanar input, e.g.
    # a prismatic skin) are deliberately KEPT: they carry no volume but they
    # do seal the boundary, and removing them punches holes in it.
    a, b, c, d = (points[tetrahedra[:, i]] for i in range(4))
    negative = numpy.einsum("ij,ij->i", b - a, numpy.cross(c - a, d - a)) < 0.0
    tetrahedra[negative, 2], tetrahedra[negative, 3] = (
        tetrahedra[negative, 3], tetrahedra[negative, 2].copy())
    return tetrahedra


def _FillWithTetgen(points, triangles, preserve_boundary, quality):
    """Constrained tetrahedralization through TetGen (pip "tetgen").

    Boundary recovery is exact. With preserve_boundary the input facets are
    kept verbatim (TetGen's -Y: Steiner points go only in the interior -
    which is all the Schoenhardt class needs, so those solids fill here);
    without it TetGen may also split boundary facets, still conforming to
    the input surface. Input vertices survive bit-identically in the leading
    rows either way; Steiner points append after them.

    TetGen itself is AGPL-licensed, which is why this backend is an explicit
    opt-in method and never part of "auto".
    """
    tetgen = _TryImportTetgen()

    generator = tetgen.TetGen(
        numpy.ascontiguousarray(points), numpy.ascontiguousarray(triangles))
    result = generator.tetrahedralize(
        order=1, nobisect=bool(preserve_boundary), quality=bool(quality))
    nodes = numpy.ascontiguousarray(result[0], dtype=numpy.float64)
    tetrahedra = numpy.asarray(result[1], dtype=numpy.int64)

    # same orientation contract as the Delaunay route: positive Jacobians
    a, b, c, d = (nodes[tetrahedra[:, i]] for i in range(4))
    negative = numpy.einsum("ij,ij->i", b - a, numpy.cross(c - a, d - a)) < 0.0
    tetrahedra[negative, 2], tetrahedra[negative, 3] = (
        tetrahedra[negative, 3], tetrahedra[negative, 2].copy())
    return nodes, tetrahedra


def FillSurfaceWithTetrahedra(surface, settings: Kratos.Parameters = None):
    """Fills a watertight 3D triangle surface with tetrahedra.

    Upstream's fill_interior is 2D-only in physicsnemo 2.2 (n = 3 raises
    NotImplementedError - "exact 3D boundary recovery is planned"), so "auto"
    tries it first and otherwise carves a Delaunay tetrahedralization by the
    winding-number sign. "tetgen" opts into a constrained tetrahedralization
    through the optional tetgen package instead - the exact-boundary-recovery
    route (see below).

    What the Delaunay fallback guarantees: every input vertex survives
    bit-identically in the leading rows, the filled volume and the boundary
    area match the input surface, and the boundary is edge-manifold - all
    three are checked, and by default a mismatch raises.

    What it does NOT guarantee: individual input facets are not preserved.
    Delaunay retriangulates planar faces with its own diagonals, so the
    boundary covers the same surface while its triangles may differ (a cube
    keeps 8 of its 12 facets). Solids that need Steiner points for boundary
    recovery - the Schoenhardt class - cannot be filled this way at all, and
    fail the validation rather than returning something wrong.

    "method": "tetgen" closes both gaps: with "preserve_boundary" (the
    default) the input facets survive VERBATIM (Steiner points are inserted
    only in the interior, which is exactly what fills a Schoenhardt-class
    solid), and the facets_preserved diagnostic asserts it. It is an
    explicit opt-in - never chosen by "auto" - both because an installed
    optional dependency must not silently change results and because TetGen
    is AGPL-licensed. Its Steiner points append after the input vertices, so
    the bit-identical-leading-rows guarantee holds there too.

    The tolerances are relative and deliberately not tighter than 1e-3: on a
    curved surface, retriangulating a non-planar boundary patch along the
    other diagonal genuinely changes the enclosed volume, so an exact match
    is only available for planar-faced input. A Schoenhardt-class failure
    misses by far more than this.

    Coplanar input (a prismatic skin, say) can leave flat slivers among the
    tetrahedra. They are kept, because they seal the boundary and dropping
    them would open it, but they are counted in the diagnostics and warned
    about - they are zero-Jacobian elements, so adapt before solving.

    The Delaunay route inserts no Steiner points: cell size follows the
    input's vertex density. For quality or size control, use "tetgen" (its
    quality pass is on by default) or adapt afterwards with
    adaptive_remeshing.RunMmgAdaptation.

    Args:
        surface: physicsnemo Mesh of TRIANGLES (cells of shape (C, 3))
            forming a closed, consistently wound surface.
        settings: Optional Kratos Parameters:
            {
                "method"                  : "auto",  // "upstream" | "delaunay" | "tetgen"
                "strict"                  : true,    // raise if validation fails
                "volume_tolerance"        : 1e-3,    // relative
                "boundary_area_tolerance" : 1e-3,    // relative
                "preserve_boundary"       : true,    // tetgen only: keep facets verbatim
                "quality"                 : true,    // tetgen only: TetGen's quality pass
                "full_output"             : false
            }

    Returns:
        A physicsnemo Mesh of tetrahedra, or (mesh, diagnostics) when
        "full_output" is set.
    """
    torch = _TryImportTorch()

    defaults = Kratos.Parameters("""{
        "method"                  : "auto",
        "strict"                  : true,
        "volume_tolerance"        : 1e-3,
        "boundary_area_tolerance" : 1e-3,
        "preserve_boundary"       : true,
        "quality"                 : true,
        "full_output"             : false
    }""")
    settings = Kratos.Parameters("{}") if settings is None else settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)

    method = settings["method"].GetString()
    if method not in ("auto", "upstream", "delaunay", "tetgen"):
        raise ValueError(
            f"Unsupported \"method\" \"{method}\". Supported: auto, upstream, delaunay, tetgen.")
    if int(surface.cells.shape[1]) != 3 or int(surface.points.shape[1]) != 3:
        raise ValueError(
            "FillSurfaceWithTetrahedra needs a 3D surface of TRIANGLES (cells of shape "
            f"(C, 3), points of shape (N, 3)); got cells {tuple(surface.cells.shape)} and "
            f"points {tuple(surface.points.shape)}. Use FillBoundaryLoop for 2D edge loops.")

    if method in ("auto", "upstream"):
        try:
            return _WithDiagnostics(_TryImportFillInterior()(surface), None, settings)
        except NotImplementedError as e:
            if method == "upstream":
                raise ValueError(
                    "\"method\": \"upstream\" needs fill_interior's 3D case, which "
                    f"physicsnemo does not implement yet ({e}). Use \"auto\" or "
                    "\"delaunay\" for the Delaunay fallback.") from e

    points = numpy.ascontiguousarray(
        numpy.asarray(surface.points.detach().cpu(), dtype=numpy.float64))
    triangles = numpy.asarray(surface.cells.detach().cpu(), dtype=numpy.int64)
    if method == "tetgen":
        filled_points, tetrahedra = _FillWithTetgen(
            points, triangles,
            settings["preserve_boundary"].GetBool(), settings["quality"].GetBool())
    else:
        filled_points, tetrahedra = points, _CarveDelaunay(surface, points)

    reference_volume, reference_area = _SurfaceMetrics(points, triangles)
    a, b, c, d = (filled_points[tetrahedra[:, i]] for i in range(4))
    filled_volume = float(
        numpy.einsum("ij,ij->i", b - a, numpy.cross(c - a, d - a)).sum() / 6.0)
    degenerate = int((numpy.abs(
        numpy.einsum("ij,ij->i", b - a, numpy.cross(c - a, d - a))
        / 6.0) <= 1e-12 * max(abs(reference_volume), 1e-300)).sum())
    boundary_faces = _BoundaryFaces(tetrahedra)
    _, filled_area = _SurfaceMetrics(filled_points, boundary_faces) if len(boundary_faces) else (0.0, 0.0)

    scale = abs(reference_volume) if abs(reference_volume) > 0.0 else 1.0
    area_scale = reference_area if reference_area > 0.0 else 1.0
    input_facets = {tuple(sorted(map(int, face))) for face in triangles}
    diagnostics = {
        "volume_ratio": filled_volume / scale,
        "boundary_area_ratio": filled_area / area_scale,
        "boundary_manifold": _IsEdgeManifold(boundary_faces),
        "tetrahedra": int(len(tetrahedra)),
        "unreferenced_points": int(len(filled_points) - len(numpy.unique(tetrahedra))),
        "degenerate_tetrahedra": degenerate,
        "steiner_points": int(len(filled_points) - len(points)),
        "facets_preserved": {tuple(map(int, face)) for face in boundary_faces} == input_facets,
    }

    failures = []
    if (method == "tetgen" and settings["preserve_boundary"].GetBool()
            and not diagnostics["facets_preserved"]):
        failures.append(
            "TetGen did not return the input facets verbatim despite preserve_boundary")
    if abs(diagnostics["volume_ratio"] - 1.0) > settings["volume_tolerance"].GetDouble():
        failures.append(
            f"filled volume {filled_volume:.9g} differs from the surface's enclosed volume "
            f"{abs(reference_volume):.9g}")
    if abs(diagnostics["boundary_area_ratio"] - 1.0) > settings["boundary_area_tolerance"].GetDouble():
        failures.append(
            f"boundary area {filled_area:.9g} differs from the input area {reference_area:.9g}")
    if not diagnostics["boundary_manifold"]:
        failures.append("the filled boundary is not edge-manifold")
    if degenerate:
        Kratos.Logger.PrintWarning(
            "FillSurfaceWithTetrahedra",
            f"{degenerate} of {len(tetrahedra)} tetrahedra are flat (zero volume). They seal "
            "the boundary and carry no volume, so the fill is geometrically correct, but they "
            "are zero-Jacobian elements: clean the mesh up with "
            "adaptive_remeshing.RunMmgAdaptation before solving on it. Coplanar input - a "
            "prismatic or otherwise structured skin - is what provokes them.")
    if failures:
        message = (
            "FillSurfaceWithTetrahedra could not fill this surface: " + "; ".join(failures) +
            ". The surface may be open, self-intersecting or need Steiner points for "
            "boundary recovery (the Schoenhardt class), which the Delaunay route cannot "
            "insert. Use \"method\": \"tetgen\" (pip install tetgen) for a constrained "
            "tetrahedralization with exact boundary recovery, or set \"strict\": false "
            "to take the result anyway.")
        if settings["strict"].GetBool():
            raise ValueError(message)
        Kratos.Logger.PrintWarning("FillSurfaceWithTetrahedra", message)

    from physicsnemo.mesh import Mesh
    mesh = Mesh(points=torch.as_tensor(filled_points), cells=torch.as_tensor(tetrahedra))
    return _WithDiagnostics(mesh, diagnostics, settings)


def _WithDiagnostics(mesh, diagnostics, settings):
    if not settings["full_output"].GetBool():
        return mesh
    if diagnostics is None:      # upstream filled it; it reports nothing of its own
        diagnostics = {"volume_ratio": None, "boundary_area_ratio": None,
                       "boundary_manifold": None,
                       "tetrahedra": int(mesh.cells.shape[0]), "unreferenced_points": None,
                       "steiner_points": None, "facets_preserved": None}
    return mesh, diagnostics


def FillModelPartWithTetrahedra(model: Kratos.Model, model_part_name: str,
                                surface_model_part: Kratos.ModelPart,
                                settings: Kratos.Parameters = None,
                                populate_settings: Kratos.Parameters = None) -> Kratos.ModelPart:
    """A surface model part in, a solvable volume model part out.

    Tessellates the surface part's Conditions through the mesh bridge, fills
    the result with FillSurfaceWithTetrahedra, and materializes it as real
    Kratos entities - the 3D counterpart of GenerateModelPart.
    """
    from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

    surface, _ = domain_mesh_builder.BuildMesh(
        surface_model_part, source_container="Conditions")
    filled = FillSurfaceWithTetrahedra(surface, settings)
    if isinstance(filled, tuple):
        filled = filled[0]
    return PopulateModelPartFromMesh(model, model_part_name, filled, populate_settings)


def RefitToImplicit(mesh, phi, iters: int = 3, bounds=None):
    """Snaps a mesh's boundary onto phi's zero set, keeping topology fixed.

    This is the DIFFERENTIABLE counterpart of GenerateImplicitDomain:
    gradients flow back through the deformed points to phi's parameters, so
    it composes with the shape-optimization layer for small perturbations.
    Large ones can invert cells (upstream warns).

    Pass `bounds` whenever the domain touches its box, or face vertices get
    dragged onto phi's interior zero set.
    """
    _TryImportTorch()
    _, _, refit_mesh_to_implicit = _TryImportGenerators()

    keywords = {"iters": int(iters)}
    if bounds is not None:
        low, high = (numpy.asarray(b, dtype=float).reshape(-1).tolist() for b in bounds)
        keywords["bounds"] = (low, high)
    return refit_mesh_to_implicit(mesh, phi, **keywords)


def PopulateModelPartFromMesh(model: Kratos.Model, model_part_name: str, mesh,
                              settings: Kratos.Parameters = None) -> Kratos.ModelPart:
    """Materializes a physicsnemo Mesh as real Kratos entities.

    The inverse of domain_mesh_builder.BuildMesh, and the step that makes a
    generated geometry solvable. Node ids are 1-based (Kratos rejects 0), so
    the mesh's 0-based connectivity is shifted here.

    Args:
        model: The Kratos Model that will own the part.
        model_part_name: Name to create. It must not already exist - Kratos
            will not silently replace entities with clashing ids.
        mesh: physicsnemo Mesh (or any object with .points/.cells).
        settings: Optional Kratos Parameters:
            {
                "entity_type"          : "Elements",   // or "Conditions"
                "element_name"         : "",           // "" = from the cell shape
                "properties_id"        : 1,
                "historical_variables" : [],           // added before the buffer
                "buffer_size"          : 1,
                "domain_size"          : 0             // 0 = the mesh's own width
            }

    Returns:
        The populated ModelPart.
    """
    torch = _TryImportTorch()

    defaults = Kratos.Parameters("""{
        "entity_type"          : "Elements",
        "element_name"         : "",
        "properties_id"        : 1,
        "historical_variables" : [],
        "buffer_size"          : 1,
        "domain_size"          : 0
    }""")
    settings = Kratos.Parameters("{}") if settings is None else settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)

    entity_type = settings["entity_type"].GetString()
    if entity_type not in ("Elements", "Conditions"):
        raise ValueError(
            f"Unsupported entity_type \"{entity_type}\". Use \"Elements\" or \"Conditions\".")
    if model.HasModelPart(model_part_name):
        raise ValueError(
            f"Model part \"{model_part_name}\" already exists; generated meshes need a "
            "fresh part (Kratos will not replace entities with clashing ids).")

    points = numpy.ascontiguousarray(
        numpy.asarray(torch.as_tensor(mesh.points).detach().cpu().numpy(), dtype=float))
    cells = numpy.ascontiguousarray(
        numpy.asarray(torch.as_tensor(mesh.cells).detach().cpu().numpy(), dtype=numpy.int64))
    if points.ndim != 2 or cells.ndim != 2:
        raise ValueError(
            f"mesh must carry (N, D) points and (C, k) cells, got {points.shape} "
            f"and {cells.shape}.")

    spatial_dimension = int(points.shape[1])
    domain_size = settings["domain_size"].GetInt() or spatial_dimension

    model_part = model.CreateModelPart(model_part_name)
    model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = domain_size
    for name in settings["historical_variables"].GetStringArray():
        model_part.AddNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(name))
    model_part.SetBufferSize(settings["buffer_size"].GetInt())

    for row, coordinates in enumerate(points):
        model_part.CreateNewNode(
            row + 1,
            float(coordinates[0]),
            float(coordinates[1]) if spatial_dimension > 1 else 0.0,
            float(coordinates[2]) if spatial_dimension > 2 else 0.0)

    properties = model_part.CreateNewProperties(settings["properties_id"].GetInt())
    entity_name = settings["element_name"].GetString() or _DefaultEntityName(
        int(cells.shape[1]), domain_size, entity_type)

    create = (model_part.CreateNewElement if entity_type == "Elements"
              else model_part.CreateNewCondition)
    for entity_id, connectivity in enumerate(cells, start=1):
        create(entity_name, entity_id,
               [int(index) + 1 for index in connectivity],   # 0-based -> 1-based
               properties)
    return model_part


def _DefaultEntityName(vertices_per_cell: int, domain_size: int, entity_type: str) -> str:
    """Kratos entity name for a simplex cell of the given width."""
    if entity_type == "Elements":
        table = {(2, 2): "Element2D2N", (3, 2): "Element2D3N",
                 (3, 3): "Element3D3N", (4, 3): "Element3D4N"}
    else:
        table = {(2, 2): "LineCondition2D2N", (3, 3): "SurfaceCondition3D3N"}
    name = table.get((vertices_per_cell, domain_size))
    if name is None:
        raise ValueError(
            f"No default {entity_type[:-1].lower()} name for cells with "
            f"{vertices_per_cell} vertices in {domain_size}D; pass \"element_name\" "
            "explicitly.")
    return name


def GenerateModelPart(model: Kratos.Model, model_part_name: str, phi, bounds, h,
                      settings: Kratos.Parameters = None,
                      populate_settings: Kratos.Parameters = None) -> Kratos.ModelPart:
    """Implicit geometry in, solvable Kratos model part out."""
    mesh = GenerateImplicitDomain(phi, bounds, h, settings)
    if isinstance(mesh, tuple):   # full_output
        mesh = mesh[0]
    return PopulateModelPartFromMesh(model, model_part_name, mesh, populate_settings)
