"""Exact NURBS geometry sampling - the isogeometric analogue of curved mode.

Kratos's NURBS geometries (IgaApplication meshes, or core
NurbsSurfaceGeometry3D / NurbsVolumeGeometry built directly) carry the exact
CAD geometry; nothing in the mesh bridge could see it, because IGA model
parts have no simplicial elements to tessellate. This module samples the
exact geometry on a structured parametric lattice - every point evaluated by
the C++ geometry's own GlobalCoordinates, never by an approximation - and
tessellates the lattice into the bridge's simplicial contract (triangles for
surfaces, Kuhn 6-tetrahedra for volumes, face-consistent across cells, so
the result is watertight by construction).

Fields in IGA live on CONTROL POINTS, and a control point is generally NOT
on the geometry - so a nodal gather makes no sense pointwise. The correct
gather is the isogeometric one: EvaluateNodalFieldOnLattice interpolates
control-point values through the geometry's own NURBS basis, obtained from
CreateQuadraturePointGeometries at the lattice's parametric coordinates -
the same the-C++-side-is-the-oracle policy the curved tessellation uses.
Sampled points are synthetic in the provenance sense: they interpolate on
gather and have no scatter-back (least-squares control-point recovery is a
fitting problem, deliberately out of scope).

Kratos stores NURBS knot vectors WITHOUT the first and last classical
repetition; the parametric domain is still [Knots[0], Knots[-1]], which is
what the lattice spans.

torch/physicsnemo are needed only by BuildNurbsMesh; sampling and field
evaluation are Kratos + numpy.
"""

import numpy

import KratosMultiphysics as Kratos


def _TryImportPhysicsNemoMesh():
    try:
        import torch
        import physicsnemo.mesh
        return torch, physicsnemo.mesh
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.nurbs_sampling requires torch and "
            "physicsnemo for BuildNurbsMesh, which could not be imported. Install "
            "them with e.g. 'pip install torch nvidia-physicsnemo'.") from e


def _HasNurbsInterface(geometry, local_dimension: int) -> bool:
    """Duck-typed: does this geometry offer everything sampling needs?

    An IGA ANALYSIS does not hand back raw NurbsSurfaceGeometry3D objects.
    Its modelers build Brep geometries that WRAP the NURBS surface, and the
    Brep's background part comes back through the Python binding as a plain
    Geometry. Both nonetheless expose the parametric span, the polynomial
    degrees, the control-point counts, exact GlobalCoordinates evaluation
    and CreateQuadraturePointGeometries - which is everything the lattice
    and the isogeometric gather actually use. Insisting on the concrete
    type would mean this bridge could sample geometries built by hand but
    not the ones a solve produces.
    """
    required = ("KnotsU", "PolynomialDegree", "PointsNumberInDirection",
                "GlobalCoordinates", "CreateQuadraturePointGeometries",
                "LocalSpaceDimension")
    if not all(hasattr(geometry, name) for name in required):
        return False
    if local_dimension == 2 and not hasattr(geometry, "KnotsV"):
        return False
    if local_dimension == 3 and not hasattr(geometry, "KnotsW"):
        return False
    try:
        return int(geometry.LocalSpaceDimension()) == local_dimension
    except Exception:
        return False


def _IsNurbsSurface(geometry) -> bool:
    return _HasNurbsInterface(geometry, 2)


def _IsNurbsVolume(geometry) -> bool:
    return _HasNurbsInterface(geometry, 3)


def _Span(knots):
    return float(knots[0]), float(knots[len(knots) - 1])


def _ResolveDivisions(divisions, count):
    if isinstance(divisions, int):
        divisions = [divisions] * count
    divisions = [int(n) for n in divisions]
    if len(divisions) != count or any(n < 1 for n in divisions):
        raise ValueError(
            f"\"divisions\" must be a positive int or {count} positive ints, "
            f"got {divisions}.")
    return divisions


def NurbsParametricLattice(geometry, divisions):
    """(local_coordinates, dims): the structured lattice in parameter space.

    local_coordinates is (M, 3) with the unused third (or second) entry 0;
    dims the per-direction point counts (divisions + 1 each). Points run
    w-slowest / u-fastest, the order the tessellations below index.
    """
    if _IsNurbsSurface(geometry):
        nu, nv = _ResolveDivisions(divisions, 2)
        u0, u1 = _Span(geometry.KnotsU())
        v0, v1 = _Span(geometry.KnotsV())
        u = numpy.linspace(u0, u1, nu + 1)
        v = numpy.linspace(v0, v1, nv + 1)
        vv, uu = numpy.meshgrid(v, u, indexing="ij")  # v-slow, u-fast
        locals_ = numpy.column_stack(
            [uu.ravel(), vv.ravel(), numpy.zeros(uu.size)])
        return locals_, (nu + 1, nv + 1)
    if _IsNurbsVolume(geometry):
        nu, nv, nw = _ResolveDivisions(divisions, 3)
        u0, u1 = _Span(geometry.KnotsU())
        v0, v1 = _Span(geometry.KnotsV())
        w0, w1 = _Span(geometry.KnotsW())
        u = numpy.linspace(u0, u1, nu + 1)
        v = numpy.linspace(v0, v1, nv + 1)
        w = numpy.linspace(w0, w1, nw + 1)
        ww, vv, uu = numpy.meshgrid(w, v, u, indexing="ij")  # w-slow, u-fast
        locals_ = numpy.column_stack([uu.ravel(), vv.ravel(), ww.ravel()])
        return locals_, (nu + 1, nv + 1, nw + 1)
    raise TypeError(
        f"Expected a NurbsSurfaceGeometry3D or NurbsVolumeGeometry, got "
        f"{type(geometry).__name__}. (Curves have no area/volume to tessellate.)")


def _LatticeTriangles(nu_points, nv_points):
    """Two consistent-diagonal triangles per lattice cell (CCW in (u, v))."""
    cells = []
    for j in range(nv_points - 1):
        for i in range(nu_points - 1):
            c00 = j * nu_points + i
            c10 = c00 + 1
            c01 = c00 + nu_points
            c11 = c01 + 1
            cells += [[c00, c10, c11], [c00, c11, c01]]
    return numpy.array(cells, dtype=numpy.int64)


# The Kuhn split: one tetrahedron per permutation of the unit steps u(1),
# v(2), w(4), each walking corner 0 -> +step -> +step -> corner 7. Every cell
# shares the same main diagonal, so shared faces split identically and the
# lattice tessellation is watertight by construction.
_KUHN_TETRAHEDRA = (
    (0, 1, 3, 7), (0, 1, 5, 7), (0, 2, 3, 7),
    (0, 2, 6, 7), (0, 4, 5, 7), (0, 4, 6, 7))


def _LatticeTetrahedra(nu_points, nv_points, nw_points):
    plane = nu_points * nv_points
    cells = []
    for k in range(nw_points - 1):
        for j in range(nv_points - 1):
            for i in range(nu_points - 1):
                corner = [
                    (k + ((code >> 2) & 1)) * plane
                    + (j + ((code >> 1) & 1)) * nu_points
                    + (i + (code & 1))
                    for code in range(8)
                ]
                cells += [[corner[c] for c in tet] for tet in _KUHN_TETRAHEDRA]
    return numpy.array(cells, dtype=numpy.int64)


def SampleNurbsGeometry(geometry, divisions):
    """(points, cells, local_coordinates): the exact geometry on the lattice.

    Every point is geometry.GlobalCoordinates at its lattice coordinate -
    exact NURBS evaluation by the C++ geometry itself. Cells are triangles
    (surface) or positively oriented tetrahedra (volume).

    Args:
        geometry: A core NurbsSurfaceGeometry3D or NurbsVolumeGeometry -
            built directly, or taken from an IgaApplication model part's
            geometries.
        divisions: Cells per parametric direction (int, or one int per
            direction).
    """
    local_coordinates, dims = NurbsParametricLattice(geometry, divisions)
    points = numpy.empty((len(local_coordinates), 3))
    for row, local in enumerate(local_coordinates):
        points[row] = geometry.GlobalCoordinates(Kratos.Vector(local))

    if len(dims) == 2:
        cells = _LatticeTriangles(*dims)
    else:
        cells = _LatticeTetrahedra(*dims)
        # a left-handed parameterization would flip every tet; restore the
        # positive-Jacobian contract the way the tetrahedral fill does
        a, b, c, d = (points[cells[:, i]] for i in range(4))
        negative = numpy.einsum("ij,ij->i", b - a, numpy.cross(c - a, d - a)) < 0.0
        cells[negative, 2], cells[negative, 3] = (
            cells[negative, 3], cells[negative, 2].copy())
    return points, cells, local_coordinates


def BuildNurbsMesh(geometry, divisions, field_specs=()):
    """A physicsnemo.mesh.Mesh of the exactly sampled geometry.

    field_specs: iterable of (variable, data_location) pairs, gathered as
    point_data through the NURBS basis (EvaluateNodalFieldOnLattice). Only
    nodal locations exist here - control points are nodes; there are no
    elements to carry element data.
    """
    torch, physicsnemo_mesh = _TryImportPhysicsNemoMesh()

    points, cells, local_coordinates = SampleNurbsGeometry(geometry, divisions)
    point_data = {}
    for variable, data_location in field_specs:
        if not str(data_location).startswith("node_"):
            raise ValueError(
                f"NURBS sampling gathers control-point (nodal) fields only; "
                f"\"{variable.Name()}\" has data_location \"{data_location}\".")
        point_data[variable.Name()] = torch.as_tensor(
            EvaluateNodalFieldOnLattice(geometry, local_coordinates, variable))
    return physicsnemo_mesh.Mesh(
        points=torch.as_tensor(points), cells=torch.as_tensor(cells),
        point_data=point_data)


def EvaluateNodalFieldOnLattice(geometry, local_coordinates, variable, step: int = 0):
    """Control-point values interpolated at the lattice points - the
    isogeometric gather.

    The basis values come from the geometry's own
    CreateQuadraturePointGeometries at the given parametric coordinates
    (each quadrature-point geometry holds exactly the nonzero control
    points and their shape function values there), so the interpolation is
    the geometry's - not a reimplementation that could drift from it.

    Returns:
        (M,) for scalar variables, (M, width) for vector ones.
    """
    quadrature_geometries = Kratos.GeometriesVector()
    integration_points = [
        [float(local[0]), float(local[1]), float(local[2]), 1.0]
        for local in local_coordinates
    ]
    geometry.CreateQuadraturePointGeometries(quadrature_geometries, 1, integration_points)
    if len(quadrature_geometries) != len(integration_points):
        raise RuntimeError(
            f"CreateQuadraturePointGeometries returned {len(quadrature_geometries)} "
            f"geometries for {len(integration_points)} lattice points; the "
            "parametric coordinates may lie outside the knot span.")

    rows = []
    for quadrature_geometry in quadrature_geometries:
        basis = numpy.asarray(quadrature_geometry.ShapeFunctionsValues())[0]
        value = None
        for weight, node in zip(basis, quadrature_geometry):
            contribution = weight * numpy.asarray(
                node.GetSolutionStepValue(variable, step), dtype=numpy.float64)
            value = contribution if value is None else value + contribution
        rows.append(value)
    return numpy.array(rows)


def SampleModelPartNurbsGeometries(model_part: Kratos.ModelPart, divisions,
                                   field_specs=()):
    """Every NURBS surface/volume geometry of a model part, sampled.

    IgaApplication modelers (NurbsGeometryModeler, CadIoModeler) leave their
    geometries in ModelPart.Geometries; this walks them and samples each
    NURBS surface or volume. Curves and non-NURBS geometries are skipped.

    Returns:
        {geometry Id: Mesh} - empty if the part carries no NURBS geometry,
        which is reported as an error since the caller asked for sampling.
    """
    meshes = {}
    for geometry in model_part.Geometries:
        if _IsNurbsSurface(geometry) or _IsNurbsVolume(geometry):
            meshes[geometry.Id] = BuildNurbsMesh(geometry, divisions, field_specs)
    if not meshes:
        raise ValueError(
            f"\"{model_part.FullName()}\" carries no NurbsSurfaceGeometry3D or "
            "NurbsVolumeGeometry (IGA modelers store them in ModelPart.Geometries).")
    return meshes
