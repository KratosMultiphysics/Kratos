"""Curved (isoparametric) tessellation of higher-order Kratos geometries.

The "subdivide" higher-order mode is straight-edged: sub-cell vertices are
the real mid-side nodes, curvature between them is lost. This module
implements ``higher_order_mode="curved"``: the parameter domain of each
quadratic geometry (Triangle6, Quadrilateral8/9, Tetrahedra10, Hexahedra27)
is refined into a dyadic lattice (``refinement_levels`` = k, ``2^k`` cells
per parametric axis) and every lattice vertex is mapped to physical space
through the exact quadratic shape functions - SYNTHETIC points that resolve
the curvature. Lattice vertices coinciding with real nodes (an exact
integer-lattice property, no epsilon) keep their node identity.

Watertightness across neighbouring curved parents is guaranteed by exact
integer *classification keys*: every lattice vertex is keyed by what it
lies on - ``(0, node_id)`` for real nodes, ``(1, A, B, i)`` for a point on
the edge between corner nodes A < B at lattice index i (measured from A),
``(2, ...)`` for points on a triangular face (sorted corner ids + permuted
integer barycentrics), ``(4, ...)`` for points on a quadrilateral face
(corner-id cycle rotated to a canonical frame + integer face coordinates in
that frame), and ``(5, entity_id, ...)`` for parent-interior points. Conforming neighbours
share their interface corner NODES, so both sides produce identical keys
for coinciding lattice points: the keys merge duplicate points into single
rows AND act as pseudo-ids for the smallest-id diagonal rules (Python
tuples order totally), so shared quadrilateral sub-faces are triangulated
identically from both sides - every interior triangle appears exactly
twice, combinatorially. No floating-point comparison is involved anywhere.

At refinement level 1 the lattice vertices of Triangle6 / Quadrilateral9 /
Hexahedra27 are exactly the real nodes and the classification keys reduce
to node ids, reproducing ``subdivide`` + ``smallest_id_diagonal``
identically (Quadrilateral8 gains one synthetic point - the parametric
center, enabling the regular 4-sub-quad pattern its missing node forbids;
curved Tetrahedra10 interior octahedra use a fixed parametric diagonal,
which never touches a parent boundary face). Serendipity solids without
the required tensor structure (Hexahedra20, Prism15, Pyramid13) and all
linear geometries fall back to the standard per-entity tessellation and
merge into the same point table - mixed meshes work, and curved mode is a
no-op for meshes without quadratic curvature.

Pure Python + numpy: this module never imports torch or physicsnemo.
"""

import functools

import numpy

import KratosMultiphysics as Kratos

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import tessellation

_GEOMETRY_TYPE = Kratos.GeometryData.KratosGeometryType

_TRIANGLE6_TYPES = (_GEOMETRY_TYPE.Kratos_Triangle2D6, _GEOMETRY_TYPE.Kratos_Triangle3D6)
_QUAD8_TYPES = (_GEOMETRY_TYPE.Kratos_Quadrilateral2D8, _GEOMETRY_TYPE.Kratos_Quadrilateral3D8)
_QUAD9_TYPES = (_GEOMETRY_TYPE.Kratos_Quadrilateral2D9, _GEOMETRY_TYPE.Kratos_Quadrilateral3D9)

CURVED_GEOMETRY_TYPES = frozenset(
    _TRIANGLE6_TYPES + _QUAD8_TYPES + _QUAD9_TYPES
    + (_GEOMETRY_TYPE.Kratos_Tetrahedra3D10, _GEOMETRY_TYPE.Kratos_Hexahedra3D27))


# --- Shape functions (transcribed from the Kratos geometry headers) ---------


def _ShapeFunctionsTriangle6(local):
    x, y = local[:, 0], local[:, 1]
    l0 = 1.0 - x - y
    return numpy.stack([
        l0 * (2.0 * l0 - 1.0), x * (2.0 * x - 1.0), y * (2.0 * y - 1.0),
        4.0 * l0 * x, 4.0 * x * y, 4.0 * y * l0], axis=1)


def _ShapeFunctionsQuadrilateral8(local):
    xi, eta = local[:, 0], local[:, 1]
    corner_xi = (-1.0, 1.0, 1.0, -1.0)
    corner_eta = (-1.0, -1.0, 1.0, 1.0)
    values = []
    for i in range(4):  # serendipity corners
        values.append(0.25 * (1.0 + xi * corner_xi[i]) * (1.0 + eta * corner_eta[i])
                      * (xi * corner_xi[i] + eta * corner_eta[i] - 1.0))
    values.append(0.5 * (1.0 - xi * xi) * (1.0 - eta))   # 4: mid (0,-1)
    values.append(0.5 * (1.0 + xi) * (1.0 - eta * eta))  # 5: mid (1,0)
    values.append(0.5 * (1.0 - xi * xi) * (1.0 + eta))   # 6: mid (0,1)
    values.append(0.5 * (1.0 - xi) * (1.0 - eta * eta))  # 7: mid (-1,0)
    return numpy.stack(values, axis=1)


def _TensorFactors(coordinate):
    return (0.5 * (coordinate - 1.0) * coordinate,   # -1 node
            0.5 * (coordinate + 1.0) * coordinate,   # +1 node
            1.0 - coordinate * coordinate)           # 0 node


def _ShapeFunctionsQuadrilateral9(local):
    fx = _TensorFactors(local[:, 0])
    fy = _TensorFactors(local[:, 1])
    order = ((0, 0), (1, 0), (1, 1), (0, 1),          # corners 0-3
             (2, 0), (1, 2), (2, 1), (0, 2), (2, 2))  # mids 4-7, center 8
    return numpy.stack([fx[i] * fy[j] for i, j in order], axis=1)


def _ShapeFunctionsTetrahedra10(local):
    x, y, z = local[:, 0], local[:, 1], local[:, 2]
    l0 = 1.0 - x - y - z
    return numpy.stack([
        l0 * (2.0 * l0 - 1.0), x * (2.0 * x - 1.0),
        y * (2.0 * y - 1.0), z * (2.0 * z - 1.0),
        4.0 * l0 * x, 4.0 * x * y, 4.0 * l0 * y,
        4.0 * l0 * z, 4.0 * x * z, 4.0 * y * z], axis=1)


# Hex27 node -> (i, j, k) tensor-factor indices (0: -1, 1: +1, 2: 0), from
# hexahedra_3d_27.h's ShapeFunctionValue case list.
_HEX27_TENSOR_ORDER = (
    (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
    (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
    (2, 0, 0), (1, 2, 0), (2, 1, 0), (0, 2, 0),
    (0, 0, 2), (1, 0, 2), (1, 1, 2), (0, 1, 2),
    (2, 0, 1), (1, 2, 1), (2, 1, 1), (0, 2, 1),
    (2, 2, 0), (2, 0, 2), (1, 2, 2), (2, 1, 2), (0, 2, 2), (2, 2, 1), (2, 2, 2))


def _ShapeFunctionsHexahedra27(local):
    fx = _TensorFactors(local[:, 0])
    fy = _TensorFactors(local[:, 1])
    fz = _TensorFactors(local[:, 2])
    return numpy.stack([fx[i] * fy[j] * fz[k] for i, j, k in _HEX27_TENSOR_ORDER], axis=1)


_SHAPE_FUNCTIONS = {}
for _type in _TRIANGLE6_TYPES:
    _SHAPE_FUNCTIONS[_type] = (_ShapeFunctionsTriangle6, 2)
for _type in _QUAD8_TYPES:
    _SHAPE_FUNCTIONS[_type] = (_ShapeFunctionsQuadrilateral8, 2)
for _type in _QUAD9_TYPES:
    _SHAPE_FUNCTIONS[_type] = (_ShapeFunctionsQuadrilateral9, 2)
_SHAPE_FUNCTIONS[_GEOMETRY_TYPE.Kratos_Tetrahedra3D10] = (_ShapeFunctionsTetrahedra10, 3)
_SHAPE_FUNCTIONS[_GEOMETRY_TYPE.Kratos_Hexahedra3D27] = (_ShapeFunctionsHexahedra27, 3)


def EvaluateShapeFunctions(geometry_type, local_coordinates) -> numpy.ndarray:
    """Shape-function values of a curved-capable geometry at local points.

    Args:
        geometry_type: One of CURVED_GEOMETRY_TYPES.
        local_coordinates: (M, d) array of local coordinates (d = 2 or 3;
            triangles/tets use their natural [0,1] coordinates, quads/hexes
            the [-1, 1] tensor coordinates).

    Returns:
        (M, n_nodes) float64, columns in the Kratos node ordering.
    """
    if geometry_type not in _SHAPE_FUNCTIONS:
        raise RuntimeError(
            f"No curved shape functions for geometry type {geometry_type}; supported: "
            f"{sorted(str(t) for t in _SHAPE_FUNCTIONS)}.")
    function, dimension = _SHAPE_FUNCTIONS[geometry_type]
    local = numpy.asarray(local_coordinates, dtype=numpy.float64).reshape(-1, dimension)
    return function(local)


# --- Parametric lattice patterns --------------------------------------------
#
# A pattern is entity-independent: lattice vertices with symbolic
# classifications (in LOCAL corner indices), exact local coordinates, the
# shape-function weight matrix, and sub-cells (some deferred until entity
# time because their triangulation is key-driven).
#
# Symbolic classifications:
#   ("real", local_node_index)
#   ("edge", corner_a, corner_b, index_from_a)             on a corner edge
#   ("tface", (ca, cb, cc), (ba, bb, bc))                  on a triangular face
#   ("face", (c0, c1, c2, c3), u, v)                       on a quad face
#   ("interior", lattice_tuple)
# Sub-cell kinds: "tri"/"tet" (final) and "quad"/"hex" (key-split at entity
# time).

_HEX_CORNER_OFFSETS = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                       (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
# quad faces of the reference hex: corner cycle + the two free axes (u, v)
_HEX_FACES = (
    ((0, 1, 2, 3), 2, 0, 0, 1),  # z = 0:  fixed axis 2 at 0, u = x, v = y
    ((4, 5, 6, 7), 2, 1, 0, 1),  # z = n
    ((0, 1, 5, 4), 1, 0, 0, 2),  # y = 0:  u = x, v = z
    ((3, 2, 6, 7), 1, 1, 0, 2),  # y = n
    ((0, 3, 7, 4), 0, 0, 1, 2),  # x = 0:  u = y, v = z
    ((1, 2, 6, 5), 0, 1, 1, 2),  # x = n
)


class _Pattern:
    def __init__(self, lattice, classifications, local_coordinates, weights, cells):
        self.lattice = lattice                    # list of integer tuples
        self.classifications = classifications    # list of symbolic tuples
        self.local_coordinates = local_coordinates  # (L, d) float64
        self.weights = weights                    # (L, n_nodes) float64
        self.cells = cells                        # list of (kind, vertex-index tuple)


def _ClassifyTensor2D(g, n, real_map, corners=(0, 1, 2, 3)):
    if g in real_map:
        return ("real", real_map[g])
    gx, gy = g
    if gy == 0:
        return ("edge", corners[0], corners[1], gx)
    if gx == n:
        return ("edge", corners[1], corners[2], gy)
    if gy == n:
        return ("edge", corners[2], corners[3], n - gx)
    if gx == 0:
        return ("edge", corners[3], corners[0], n - gy)
    return ("face", corners, gx, gy)


@functools.lru_cache(maxsize=None)
def _QuadrilateralPattern(geometry_type, levels):
    n = 2 ** levels
    half = n // 2
    real_map = {(0, 0): 0, (n, 0): 1, (n, n): 2, (0, n): 3,
                (half, 0): 4, (n, half): 5, (half, n): 6, (0, half): 7}
    if geometry_type in _QUAD9_TYPES:
        real_map[(half, half)] = 8  # Quad8: the center stays synthetic

    lattice, index_of = [], {}
    for gx in range(n + 1):
        for gy in range(n + 1):
            index_of[(gx, gy)] = len(lattice)
            lattice.append((gx, gy))
    classifications = [_ClassifyTensor2D(g, n, real_map) for g in lattice]
    local = numpy.array([(2.0 * gx / n - 1.0, 2.0 * gy / n - 1.0) for gx, gy in lattice])
    weights = EvaluateShapeFunctions(geometry_type, local)

    cells = []
    for gx in range(n):
        for gy in range(n):
            cells.append(("quad", (index_of[(gx, gy)], index_of[(gx + 1, gy)],
                                   index_of[(gx + 1, gy + 1)], index_of[(gx, gy + 1)])))
    return _Pattern(lattice, classifications, local, weights, cells)


@functools.lru_cache(maxsize=None)
def _TrianglePattern(geometry_type, levels):
    n = 2 ** levels
    half = n // 2
    real_map = {(0, 0): 0, (n, 0): 1, (0, n): 2,
                (half, 0): 3, (half, half): 4, (0, half): 5}

    lattice, index_of = [], {}
    for a in range(n + 1):
        for b in range(n + 1 - a):
            index_of[(a, b)] = len(lattice)
            lattice.append((a, b))

    def classify(g):
        if g in real_map:
            return ("real", real_map[g])
        a, b = g
        if b == 0:
            return ("edge", 0, 1, a)
        if a + b == n:
            return ("edge", 1, 2, b)
        if a == 0:
            return ("edge", 2, 0, n - b)
        return ("tface", (0, 1, 2), (n - a - b, a, b))

    classifications = [classify(g) for g in lattice]
    local = numpy.array([(a / n, b / n) for a, b in lattice])
    weights = EvaluateShapeFunctions(geometry_type, local)

    cells = []
    for a in range(n):
        for b in range(n - a):
            cells.append(("tri", (index_of[(a, b)], index_of[(a + 1, b)], index_of[(a, b + 1)])))
            if a + b <= n - 2:
                cells.append(("tri", (index_of[(a + 1, b)], index_of[(a + 1, b + 1)],
                                      index_of[(a, b + 1)])))
    return _Pattern(lattice, classifications, local, weights, cells)


@functools.lru_cache(maxsize=None)
def _HexahedronPattern(levels):
    geometry_type = _GEOMETRY_TYPE.Kratos_Hexahedra3D27
    n = 2 ** levels
    half = n // 2
    # real nodes at every lattice position with components in {0, half, n}
    factor_position = {0: 0, 1: n, 2: half}  # tensor index -> lattice component
    real_map = {}
    for node, (i, j, k) in enumerate(_HEX27_TENSOR_ORDER):
        real_map[(factor_position[i], factor_position[j], factor_position[k])] = node

    lattice, index_of = [], {}
    for gx in range(n + 1):
        for gy in range(n + 1):
            for gz in range(n + 1):
                index_of[(gx, gy, gz)] = len(lattice)
                lattice.append((gx, gy, gz))

    hex_edges = tuple((a, b) for a in range(8) for b in range(a + 1, 8)
                      if sum(abs(_HEX_CORNER_OFFSETS[a][c] - _HEX_CORNER_OFFSETS[b][c])
                             for c in range(3)) == 1)

    def classify(g):
        if g in real_map:
            return ("real", real_map[g])
        extreme = [c for c in range(3) if g[c] in (0, n)]
        if len(extreme) == 2:  # on a corner edge
            free_axis = ({0, 1, 2} - set(extreme)).pop()
            for a, b in hex_edges:
                pa = tuple(_HEX_CORNER_OFFSETS[a][c] * n for c in range(3))
                pb = tuple(_HEX_CORNER_OFFSETS[b][c] * n for c in range(3))
                if not all(g[c] == pa[c] for c in extreme):
                    continue
                if pa[free_axis] == 0 and pb[free_axis] == n:
                    return ("edge", a, b, g[free_axis])
                if pb[free_axis] == 0 and pa[free_axis] == n:
                    return ("edge", b, a, g[free_axis])
            raise RuntimeError(f"unclassifiable edge point {g}")  # pragma: no cover
        if len(extreme) == 1:
            axis = extreme[0]
            side = 1 if g[axis] == n else 0
            for cycle, face_axis, face_side, u_axis, v_axis in _HEX_FACES:
                if face_axis == axis and face_side == side:
                    return ("face", cycle, g[u_axis], g[v_axis])
            raise RuntimeError(f"unclassifiable face point {g}")  # pragma: no cover
        return ("interior", g)

    classifications = [classify(g) for g in lattice]
    local = numpy.array([(2.0 * gx / n - 1.0, 2.0 * gy / n - 1.0, 2.0 * gz / n - 1.0)
                         for gx, gy, gz in lattice])
    weights = EvaluateShapeFunctions(geometry_type, local)

    cells = []
    for gx in range(n):
        for gy in range(n):
            for gz in range(n):
                corners = tuple(index_of[(gx + dx, gy + dy, gz + dz)]
                                for dx, dy, dz in _HEX_CORNER_OFFSETS)
                cells.append(("hex", corners))
    return _Pattern(lattice, classifications, local, weights, cells)


@functools.lru_cache(maxsize=None)
def _TetrahedronPattern(levels):
    geometry_type = _GEOMETRY_TYPE.Kratos_Tetrahedra3D10
    n = 2 ** levels
    half = n // 2
    # integer barycentric-style coordinates (x, y, z) with x + y + z <= n
    real_map = {(0, 0, 0): 0, (n, 0, 0): 1, (0, n, 0): 2, (0, 0, n): 3,
                (half, 0, 0): 4, (half, half, 0): 5, (0, half, 0): 6,
                (0, 0, half): 7, (half, 0, half): 8, (0, half, half): 9}

    lattice, index_of = [], {}

    def vertex(g):
        if g not in index_of:
            index_of[g] = len(lattice)
            lattice.append(g)
        return index_of[g]

    cells = []

    def refine(v0, v1, v2, v3, depth):
        if depth == 0:
            cells.append(("tet", (vertex(v0), vertex(v1), vertex(v2), vertex(v3))))
            return
        mid = lambda a, b: tuple((a[c] + b[c]) // 2 for c in range(3))
        m01, m12, m02 = mid(v0, v1), mid(v1, v2), mid(v0, v2)
        m03, m13, m23 = mid(v0, v3), mid(v1, v3), mid(v2, v3)
        refine(v0, m01, m02, m03, depth - 1)
        refine(m01, v1, m12, m13, depth - 1)
        refine(m02, m12, v2, m23, depth - 1)
        refine(m03, m13, m23, v3, depth - 1)
        # interior octahedron: fixed parametric diagonal m01-m23 (never on a
        # parent boundary face, so cross-parent watertightness is unaffected)
        refine(m01, m12, m02, m23, depth - 1)
        refine(m01, m12, m13, m23, depth - 1)
        refine(m01, m02, m03, m23, depth - 1)
        refine(m01, m03, m13, m23, depth - 1)

    refine((0, 0, 0), (n, 0, 0), (0, n, 0), (0, 0, n), levels)

    tet_edges = ((0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3))
    corner_positions = {0: (0, 0, 0), 1: (n, 0, 0), 2: (0, n, 0), 3: (0, 0, n)}
    tet_faces = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))

    def barycentrics(g):
        return (n - g[0] - g[1] - g[2], g[0], g[1], g[2])  # wrt corners 0..3

    def classify(g):
        if g in real_map:
            return ("real", real_map[g])
        bary = barycentrics(g)
        zero_corners = [c for c in range(4) if bary[c] == 0]
        if len(zero_corners) == 2:  # on the edge between the two NON-zero corners
            a, b = [c for c in range(4) if bary[c] > 0]
            for ea, eb in tet_edges:
                if {ea, eb} == {a, b}:
                    return ("edge", ea, eb, bary[eb])
            raise RuntimeError(f"unclassifiable tet edge point {g}")  # pragma: no cover
        if len(zero_corners) == 1:
            face = tuple(c for c in range(4) if c != zero_corners[0])
            return ("tface", face, tuple(bary[c] for c in face))
        return ("interior", g)

    classifications = [classify(g) for g in lattice]
    local = numpy.array([(g[0] / n, g[1] / n, g[2] / n) for g in lattice])
    weights = EvaluateShapeFunctions(geometry_type, local)
    return _Pattern(lattice, classifications, local, weights, cells)


def _GetPattern(geometry_type, levels):
    if geometry_type in _TRIANGLE6_TYPES:
        return _TrianglePattern(geometry_type, levels)
    if geometry_type in _QUAD8_TYPES or geometry_type in _QUAD9_TYPES:
        return _QuadrilateralPattern(geometry_type, levels)
    if geometry_type == _GEOMETRY_TYPE.Kratos_Tetrahedra3D10:
        return _TetrahedronPattern(levels)
    if geometry_type == _GEOMETRY_TYPE.Kratos_Hexahedra3D27:
        return _HexahedronPattern(levels)
    raise RuntimeError(f"Geometry type {geometry_type} has no curved pattern.")  # pragma: no cover


# --- Entity-time key construction and key-driven splitting ------------------


def _CanonicalQuadFaceKey(cycle_ids, u, v, n):
    """Canonical key of a point on a quadrilateral face lattice.

    cycle_ids: the face's 4 corner node ids in cyclic order; (u, v) the
    point's lattice coordinates in the frame with origin at cycle position
    0, u toward position 1, v toward position 3. Both parents of a shared
    face compute identical keys because the canonical frame depends only on
    the (shared) corner ids.
    """
    corner_uv = ((0, 0), (n, 0), (n, n), (0, n))
    r = min(range(4), key=lambda i: cycle_ids[i])
    direction = 1 if cycle_ids[(r + 1) % 4] < cycle_ids[(r - 1) % 4] else -1
    origin = corner_uv[r]
    e_u = tuple((corner_uv[(r + direction) % 4][c] - origin[c]) // n for c in range(2))
    e_v = tuple((corner_uv[(r - direction) % 4][c] - origin[c]) // n for c in range(2))
    du, dv = u - origin[0], v - origin[1]
    u_prime = du * e_u[0] + dv * e_u[1]
    v_prime = du * e_v[0] + dv * e_v[1]
    # category tag 4 (quad face) - distinct from 2 (triangular face) so the
    # two 7-tuple layouts can never collide
    return (4, cycle_ids[r], cycle_ids[(r + direction) % 4],
            cycle_ids[(r + 2 * direction) % 4], cycle_ids[(r - direction) % 4],
            u_prime, v_prime)


def _BuildKeys(pattern, node_ids, entity_id, n):
    """Classification key of every lattice vertex, for one entity."""
    keys = []
    for classification in pattern.classifications:
        kind = classification[0]
        if kind == "real":
            keys.append((0, node_ids[classification[1]]))
        elif kind == "edge":
            _, a, b, index = classification
            id_a, id_b = node_ids[a], node_ids[b]
            if id_a < id_b:
                keys.append((1, id_a, id_b, index))
            else:
                keys.append((1, id_b, id_a, n - index))
        elif kind == "tface":
            _, corners, barycentrics = classification
            pairs = sorted(zip((node_ids[c] for c in corners), barycentrics))
            keys.append((2,) + tuple(p[0] for p in pairs) + tuple(p[1] for p in pairs))
        elif kind == "face":
            _, cycle, u, v = classification
            keys.append(_CanonicalQuadFaceKey(tuple(node_ids[c] for c in cycle), u, v, n))
        else:  # interior (tag 5: above every shareable category)
            keys.append((5, entity_id) + tuple(classification[1]))
    return keys


def _SplitQuadByKeys(vertices, keys):
    """Splits a quad (CCW vertex indices) along the smallest-key diagonal."""
    smallest = min(range(4), key=lambda i: keys[vertices[i]])
    if smallest % 2 == 0:
        return [(vertices[0], vertices[1], vertices[2]), (vertices[0], vertices[2], vertices[3])]
    return [(vertices[0], vertices[1], vertices[3]), (vertices[1], vertices[2], vertices[3])]


def _SplitHexByKeys(vertices, keys):
    """Dompierre split of a sub-hex, driven by the vertex classification keys."""
    corner_keys = [keys[v] for v in vertices]
    index_of_key = {key: vertices[i] for i, key in enumerate(corner_keys)}
    return [tuple(index_of_key[key] for key in simplex)
            for simplex in tessellation._TessellateHexahedronSmallestId(corner_keys)]


# --- Container tessellation --------------------------------------------------


class CurvedTessellationResult:
    """Point/cell arrays of a curved container tessellation.

    Points are ordered real-first (ascending node id) then synthetic
    (deterministic key order); synthetic entries carry parent, local
    coordinates and shape-function weights for interpolation.
    """

    def __init__(self, point_node_ids, point_coordinates, simplex_cells, cell_provenance,
                 synthetic_parent_ids, synthetic_local_coordinates,
                 synthetic_node_ids, synthetic_weights):
        self.point_node_ids = point_node_ids
        self.point_coordinates = point_coordinates
        self.simplex_cells = simplex_cells
        self.cell_provenance = cell_provenance
        self.synthetic_parent_ids = synthetic_parent_ids
        self.synthetic_local_coordinates = synthetic_local_coordinates
        self.synthetic_node_ids = synthetic_node_ids
        self.synthetic_weights = synthetic_weights


def TessellateContainerCurved(container, node_coordinates, refinement_levels: int = 2):
    """Curved tessellation of every entity of a container.

    Args:
        container: Kratos Elements or Conditions container.
        node_coordinates: dict mapping node id -> numpy (3,) coordinates.
        refinement_levels: dyadic parametric refinement depth k >= 1
            (2^k cells per parametric axis).

    Returns:
        CurvedTessellationResult (see the module docstring for the
        point-ordering and watertightness contracts).
    """
    refinement_levels = int(refinement_levels)
    if refinement_levels < 1:
        raise ValueError(f"refinement_levels must be >= 1, got {refinement_levels}.")
    n = 2 ** refinement_levels

    point_index_by_key = {}
    point_records = []  # per point: (key, node_id or -1, coords, synthetic data or None)
    simplices = []
    cell_provenance = []

    def get_point(key, coordinates, node_id=-1, synthetic=None):
        index = point_index_by_key.get(key)
        if index is None:
            index = len(point_records)
            point_index_by_key[key] = index
            point_records.append((key, node_id, coordinates, synthetic))
        return index

    for entity in container:
        geometry = entity.GetGeometry()
        geometry_type = geometry.GetGeometryType()
        node_ids = [node.Id for node in geometry]
        entity_id = entity.Id

        if geometry_type not in CURVED_GEOMETRY_TYPES:
            # linear geometries and unsupported serendipity solids: the
            # standard per-entity tessellation, merged via (0, id) keys
            for sub_index, simplex in enumerate(tessellation.TessellateEntity(
                    geometry_type, node_ids, node_coordinates,
                    "smallest_id_diagonal", "reduce")):
                simplices.append(tuple(
                    get_point((0, node_id), node_coordinates[node_id], node_id)
                    for node_id in simplex))
                cell_provenance.append((entity_id, sub_index))
            continue

        pattern = _GetPattern(geometry_type, refinement_levels)
        keys = _BuildKeys(pattern, node_ids, entity_id, n)
        entity_coordinates = numpy.stack([node_coordinates[node_id] for node_id in node_ids])
        physical = pattern.weights @ entity_coordinates  # (L, 3)

        vertex_to_point = []
        for lattice_index, key in enumerate(keys):
            classification = pattern.classifications[lattice_index]
            if classification[0] == "real":
                node_id = node_ids[classification[1]]
                vertex_to_point.append(get_point(key, node_coordinates[node_id], node_id))
            else:
                synthetic = (entity_id,
                             pattern.local_coordinates[lattice_index],
                             node_ids,
                             pattern.weights[lattice_index])
                vertex_to_point.append(get_point(key, physical[lattice_index], -1, synthetic))

        point_keys = {point: point_records[point][0] for point in vertex_to_point}
        sub_index = 0
        for kind, vertices in pattern.cells:
            mapped = tuple(vertex_to_point[v] for v in vertices)
            if kind in ("tri", "tet"):
                cell_simplices = [mapped]
            elif kind == "quad":
                cell_simplices = _SplitQuadByKeys(mapped, point_keys)
            else:  # hex
                cell_simplices = _SplitHexByKeys(mapped, point_keys)
            for simplex in cell_simplices:
                simplices.append(simplex)
                cell_provenance.append((entity_id, sub_index))
                sub_index += 1

    # --- reorder: real points ascending by node id, then synthetic by key ---
    real = [(record[1], index) for index, record in enumerate(point_records) if record[1] >= 0]
    synthetic = [(record[0], index) for index, record in enumerate(point_records) if record[1] < 0]
    real.sort()
    synthetic.sort()
    new_order = [index for _, index in real] + [index for _, index in synthetic]
    remap = numpy.empty(len(point_records), dtype=numpy.int64)
    remap[new_order] = numpy.arange(len(point_records))

    point_node_ids = numpy.array(
        [point_records[index][1] for index in new_order], dtype=numpy.int64)
    point_coordinates = numpy.array(
        [point_records[index][2] for index in new_order], dtype=numpy.float64)
    simplex_cells = remap[numpy.asarray(simplices, dtype=numpy.int64)]

    synthetic_records = [point_records[index][3] for _, index in synthetic]
    if synthetic_records:
        n_max = max(len(record[2]) for record in synthetic_records)
        parent_ids = numpy.array([record[0] for record in synthetic_records], dtype=numpy.int64)
        local_coordinates = numpy.zeros((len(synthetic_records), 3))
        node_id_rows = numpy.empty((len(synthetic_records), n_max), dtype=numpy.int64)
        weight_rows = numpy.zeros((len(synthetic_records), n_max))
        for row, (_, local, ids, weights) in enumerate(synthetic_records):
            local_coordinates[row, :len(local)] = local
            node_id_rows[row, :len(ids)] = ids
            node_id_rows[row, len(ids):] = ids[0]  # padding: repeated id, weight 0
            weight_rows[row, :len(weights)] = weights
    else:
        parent_ids = local_coordinates = node_id_rows = weight_rows = None

    return CurvedTessellationResult(
        point_node_ids=point_node_ids,
        point_coordinates=point_coordinates,
        simplex_cells=simplex_cells,
        cell_provenance=numpy.array(cell_provenance, dtype=numpy.int64).reshape(-1, 2),
        synthetic_parent_ids=parent_ids,
        synthetic_local_coordinates=local_coordinates,
        synthetic_node_ids=node_id_rows,
        synthetic_weights=weight_rows)
