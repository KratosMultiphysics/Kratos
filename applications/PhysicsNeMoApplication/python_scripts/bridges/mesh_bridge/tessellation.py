"""Tessellation of Kratos geometries into simplices.

physicsnemo.mesh handles strictly simplicial meshes (triangles/tetrahedra), so
non-simplex Kratos geometries have to be decomposed. All decompositions are
built from existing Kratos nodes only: no synthetic points are ever
introduced, every simplex vertex is an original node of the entity.

Two tessellation modes are supported (``mode``):

- ``"smallest_id_diagonal"`` (default): every quadrilateral face is split
  along the diagonal that passes through the face's smallest global node id
  (Dompierre, Labbe, Vallet & Camarero, "How to Subdivide Pyramids, Prisms
  and Hexahedra into Tetrahedra", IMR 1999). The diagonal choice depends
  only on the global ids of the face nodes, so two entities sharing a face
  always triangulate it identically: the tessellation of any conforming mesh
  (structured or not, and including quadrilateral surface conditions on
  hexahedron faces) is watertight. Hexahedra decompose into 5 or 6
  tetrahedra depending on the resulting diagonal configuration.

- ``"fan"``: the legacy fixed tables (6-tetrahedra fan around the
  hexahedron 0-6 diagonal, shortest-diagonal quadrilateral split). Face
  triangulations are only consistent between neighbours whose local node
  numbering is translationally consistent (e.g. structured grids).

Higher-order handling (``higher_order_mode``):

- ``"reduce"`` (default): higher-order geometries are reduced to their
  linear corner sub-geometry (Kratos orders corner nodes first); mid-side
  and interior nodes - and their field values - are dropped.

- ``"subdivide"``: geometries with a full set of mid-side/interior nodes are
  subdivided into linear sub-entities through those real nodes: Triangle6 ->
  4 triangles, Quadrilateral8 -> 6 triangles, Quadrilateral9 -> 8 triangles,
  Tetrahedra10 -> 8 tetrahedra, Hexahedra27 -> 8 sub-hexahedra tessellated
  per the active mode. Serendipity types lacking the interior nodes needed
  for a conforming subdivision (Hexahedra20, Prism15, Pyramid13) fall back
  to corner reduction with a one-time warning. The subdivision is a
  straight-edged approximation of the curved geometry: sub-entity vertices
  interpolate the true geometry, curvature between them is lost.

Pure Python + numpy: this module never imports torch or physicsnemo.
"""

import itertools

import numpy

import KratosMultiphysics as Kratos

_GEOMETRY_TYPE = Kratos.GeometryData.KratosGeometryType

_TESSELLATION_MODES = ("smallest_id_diagonal", "fan")
_HIGHER_ORDER_MODES = ("reduce", "subdivide")

# Number of corner nodes per supported geometry type (corners come first in
# the Kratos node ordering, so higher-order types just truncate).
_CORNER_COUNT = {
    _GEOMETRY_TYPE.Kratos_Triangle2D3: 3,
    _GEOMETRY_TYPE.Kratos_Triangle2D6: 3,
    _GEOMETRY_TYPE.Kratos_Triangle3D3: 3,
    _GEOMETRY_TYPE.Kratos_Triangle3D6: 3,
    _GEOMETRY_TYPE.Kratos_Quadrilateral2D4: 4,
    _GEOMETRY_TYPE.Kratos_Quadrilateral2D8: 4,
    _GEOMETRY_TYPE.Kratos_Quadrilateral2D9: 4,
    _GEOMETRY_TYPE.Kratos_Quadrilateral3D4: 4,
    _GEOMETRY_TYPE.Kratos_Quadrilateral3D8: 4,
    _GEOMETRY_TYPE.Kratos_Quadrilateral3D9: 4,
    _GEOMETRY_TYPE.Kratos_Tetrahedra3D4: 4,
    _GEOMETRY_TYPE.Kratos_Tetrahedra3D10: 4,
    _GEOMETRY_TYPE.Kratos_Hexahedra3D8: 8,
    _GEOMETRY_TYPE.Kratos_Hexahedra3D20: 8,
    _GEOMETRY_TYPE.Kratos_Hexahedra3D27: 8,
    _GEOMETRY_TYPE.Kratos_Prism3D6: 6,
    _GEOMETRY_TYPE.Kratos_Prism3D15: 6,
    _GEOMETRY_TYPE.Kratos_Pyramid3D5: 5,
    _GEOMETRY_TYPE.Kratos_Pyramid3D13: 5,
}

_TRIANGLE_TYPES = {
    _GEOMETRY_TYPE.Kratos_Triangle2D3, _GEOMETRY_TYPE.Kratos_Triangle2D6,
    _GEOMETRY_TYPE.Kratos_Triangle3D3, _GEOMETRY_TYPE.Kratos_Triangle3D6,
}
_TETRAHEDRON_TYPES = {_GEOMETRY_TYPE.Kratos_Tetrahedra3D4, _GEOMETRY_TYPE.Kratos_Tetrahedra3D10}
_QUADRILATERAL_TYPES = {
    _GEOMETRY_TYPE.Kratos_Quadrilateral2D4, _GEOMETRY_TYPE.Kratos_Quadrilateral2D8,
    _GEOMETRY_TYPE.Kratos_Quadrilateral2D9, _GEOMETRY_TYPE.Kratos_Quadrilateral3D4,
    _GEOMETRY_TYPE.Kratos_Quadrilateral3D8, _GEOMETRY_TYPE.Kratos_Quadrilateral3D9,
}
_HEXAHEDRON_TYPES = {
    _GEOMETRY_TYPE.Kratos_Hexahedra3D8, _GEOMETRY_TYPE.Kratos_Hexahedra3D20,
    _GEOMETRY_TYPE.Kratos_Hexahedra3D27,
}
_PRISM_TYPES = {_GEOMETRY_TYPE.Kratos_Prism3D6, _GEOMETRY_TYPE.Kratos_Prism3D15}
_PYRAMID_TYPES = {_GEOMETRY_TYPE.Kratos_Pyramid3D5, _GEOMETRY_TYPE.Kratos_Pyramid3D13}

# --- Legacy "fan" tables ----------------------------------------------------

# 6-tetrahedra fan around the 0-6 diagonal (Kratos hex ordering: 0-3 bottom
# counter-clockwise, 4-7 top counter-clockwise above them).
_HEXAHEDRON_FAN_TABLE = (
    (0, 1, 2, 6),
    (0, 2, 3, 6),
    (0, 3, 7, 6),
    (0, 7, 4, 6),
    (0, 4, 5, 6),
    (0, 5, 1, 6),
)

# Prism/wedge (0-2 bottom triangle, 3-5 top triangle).
_PRISM_FAN_TABLE = (
    (0, 1, 2, 3),
    (1, 2, 3, 4),
    (2, 3, 4, 5),
)

# Pyramid (0-3 base quadrilateral, 4 apex), split along base diagonal 0-2.
_PYRAMID_FAN_TABLE = (
    (0, 1, 2, 4),
    (0, 2, 3, 4),
)

# --- Smallest-id (Dompierre) machinery --------------------------------------


def _BuildHexahedronRotations():
    """All 24 orientation-preserving relabellings of a hexahedron, grouped by
    the local corner index each relabelling brings to position 0.

    A permutation p reorders the corners as (corners[p[0]], ..., corners[p[7]])
    while keeping valid Kratos hexahedron connectivity and orientation. There
    are exactly 3 per leading corner (the rotations about its main diagonal).
    """
    reference = numpy.array(
        [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
         [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]])
    centered = reference - 0.5
    by_leading_corner = {}
    for axes in itertools.permutations(range(3)):
        for signs in itertools.product((1.0, -1.0), repeat=3):
            rotation = numpy.zeros((3, 3))
            for row, (axis, sign) in enumerate(zip(axes, signs)):
                rotation[row, axis] = sign
            if numpy.linalg.det(rotation) < 0.0:
                continue
            rotated = centered @ rotation.T + 0.5
            destination = [int(numpy.abs(reference - point).sum(axis=1).argmin()) for point in rotated]
            permutation = [0] * 8
            for old_index, new_index in enumerate(destination):
                permutation[new_index] = old_index
            by_leading_corner.setdefault(permutation[0], []).append(tuple(permutation))
    return by_leading_corner


_HEXAHEDRON_ROTATIONS = _BuildHexahedronRotations()

# Orientation-preserving prism relabellings bringing each corner to local 0
# (vertical-axis rotations and 180-degree horizontal flips).
_PRISM_ROTATIONS = {
    0: (0, 1, 2, 3, 4, 5),
    1: (1, 2, 0, 4, 5, 3),
    2: (2, 0, 1, 5, 3, 4),
    3: (3, 5, 4, 0, 2, 1),
    4: (4, 3, 5, 1, 0, 2),
    5: (5, 4, 3, 2, 1, 0),
}

# The three hexahedron faces incident to local corner 6, the corner opposite
# the (rotated-to-0) smallest node id. Faces containing 0 always take their
# diagonal through 0, the global minimum, so only these three faces vary.
_HEXAHEDRON_FACES_AT_SIX = ((1, 2, 6, 5), (2, 3, 7, 6), (4, 5, 6, 7))

# Dompierre decomposition tables after rotating the smallest global node id
# to local 0. Keyed by whether the smallest-id diagonal of each face in
# _HEXAHEDRON_FACES_AT_SIX passes through corner 6. Only 4 of the 8 keys are
# canonical; the others are reached from them by a rotation about the 0-6
# diagonal, which _TessellateHexahedronSmallestId searches over.
_HEXAHEDRON_SMALLEST_ID_TABLES = {
    (False, False, False): (  # 5 tetrahedra: corner tets at 1, 3, 4, 6 + core
        (0, 1, 2, 5), (0, 2, 3, 7), (0, 4, 5, 7), (2, 5, 6, 7), (0, 2, 7, 5)),
    (True, False, False): (
        (0, 1, 2, 6), (0, 1, 6, 5), (0, 5, 7, 4), (0, 5, 6, 7), (0, 2, 7, 6), (0, 2, 3, 7)),
    (True, True, False): (
        (0, 1, 2, 6), (0, 1, 6, 5), (0, 2, 3, 6), (0, 3, 7, 6), (0, 5, 6, 7), (0, 5, 7, 4)),
    (True, True, True): _HEXAHEDRON_FAN_TABLE,
}

# Quad face of the prism not containing local corner 0, after rotation. Its
# diagonal is the only free choice (both faces at 0 go through 0).
_PRISM_FACE_AT_ZERO_OPPOSITE = (1, 2, 5, 4)
_PRISM_SMALLEST_ID_TABLES = {
    True: ((0, 1, 2, 5), (0, 1, 5, 4), (0, 4, 5, 3)),   # diagonal 1-5
    False: ((0, 1, 2, 4), (0, 2, 5, 4), (0, 4, 5, 3)),  # diagonal 2-4
}


def _DiagonalThroughCorner(face_ids, corner_position):
    """Whether the smallest-id diagonal of quad face (a,b,c,d) - diagonals
    a-c (positions 0,2) and b-d (positions 1,3) - contains corner_position."""
    smallest = min(range(4), key=lambda k: face_ids[k])
    return smallest % 2 == corner_position % 2


def _TessellateHexahedronSmallestId(corners):
    smallest = min(range(8), key=lambda i: corners[i])
    for permutation in _HEXAHEDRON_ROTATIONS[smallest]:
        rotated = tuple(corners[i] for i in permutation)
        flags = tuple(
            _DiagonalThroughCorner(tuple(rotated[i] for i in face), face.index(6))
            for face in _HEXAHEDRON_FACES_AT_SIX)
        table = _HEXAHEDRON_SMALLEST_ID_TABLES.get(flags)
        if table is not None:
            return [tuple(rotated[i] for i in simplex) for simplex in table]
    raise RuntimeError(  # pragma: no cover - one of the 3 rotations always matches
        f"No canonical diagonal configuration found for hexahedron corners {corners}.")


def _TessellatePrismSmallestId(corners):
    smallest = min(range(6), key=lambda i: corners[i])
    rotated = tuple(corners[i] for i in _PRISM_ROTATIONS[smallest])
    face_ids = tuple(rotated[i] for i in _PRISM_FACE_AT_ZERO_OPPOSITE)
    table = _PRISM_SMALLEST_ID_TABLES[_DiagonalThroughCorner(face_ids, 0)]
    return [tuple(rotated[i] for i in simplex) for simplex in table]


def _TessellatePyramidSmallestId(corners):
    if min(range(4), key=lambda i: corners[i]) % 2 == 0:
        table = ((0, 1, 2, 4), (0, 2, 3, 4))  # base diagonal 0-2
    else:
        table = ((0, 1, 3, 4), (1, 2, 3, 4))  # base diagonal 1-3
    return [tuple(corners[i] for i in simplex) for simplex in table]


def _SquaredDistance(coordinates, node_id_a, node_id_b):
    d = coordinates[node_id_a] - coordinates[node_id_b]
    return float(d @ d)


def _SplitQuadrilateral(corners, node_coordinates, mode):
    if mode == "smallest_id_diagonal":
        through_zero_two = min(range(4), key=lambda i: corners[i]) % 2 == 0
    else:  # fan: split along the shortest diagonal to avoid slivers
        through_zero_two = _SquaredDistance(node_coordinates, corners[0], corners[2]) <= \
            _SquaredDistance(node_coordinates, corners[1], corners[3])
    if through_zero_two:
        return [(corners[0], corners[1], corners[2]),
                (corners[0], corners[2], corners[3])]
    return [(corners[0], corners[1], corners[3]),
            (corners[1], corners[2], corners[3])]


# --- Higher-order subdivision tables (Kratos node orderings) ----------------

_TRIANGLE6_TYPES = {_GEOMETRY_TYPE.Kratos_Triangle2D6, _GEOMETRY_TYPE.Kratos_Triangle3D6}
_QUAD8_TYPES = {_GEOMETRY_TYPE.Kratos_Quadrilateral2D8, _GEOMETRY_TYPE.Kratos_Quadrilateral3D8}
_QUAD9_TYPES = {_GEOMETRY_TYPE.Kratos_Quadrilateral2D9, _GEOMETRY_TYPE.Kratos_Quadrilateral3D9}
_UNSUBDIVIDABLE_TYPES = {
    _GEOMETRY_TYPE.Kratos_Hexahedra3D20,
    _GEOMETRY_TYPE.Kratos_Prism3D15,
    _GEOMETRY_TYPE.Kratos_Pyramid3D13,
}

# Triangle6: mid-edge nodes 3 (0-1), 4 (1-2), 5 (2-0).
_TRIANGLE6_SUBDIVISION = ((0, 3, 5), (3, 1, 4), (5, 4, 2), (3, 4, 5))

# Quadrilateral8/9: mid-edge nodes 4 (0-1), 5 (1-2), 6 (2-3), 7 (3-0); 8 center.
_QUAD9_SUB_QUADS = ((0, 4, 8, 7), (4, 1, 5, 8), (8, 5, 2, 6), (7, 8, 6, 3))
_QUAD8_CORNER_TRIANGLES = ((0, 4, 7), (1, 5, 4), (2, 6, 5), (3, 7, 6))
_QUAD8_INNER_QUAD = (4, 5, 6, 7)

# Tetrahedra10: mid-edge nodes 4 (0-1), 5 (1-2), 6 (0-2), 7 (0-3), 8 (1-3),
# 9 (2-3). Four corner tets + the central octahedron, split around one of its
# three diagonals (opposite mid-node pairs).
_TET10_CORNER_TETS = ((0, 4, 6, 7), (4, 1, 5, 8), (6, 5, 2, 9), (7, 8, 9, 3))
_TET10_OCTAHEDRON_DIAGONALS = ((4, 9), (5, 7), (6, 8))
_TET10_OCTAHEDRON_TABLES = {
    (4, 9): ((4, 5, 6, 9), (4, 6, 7, 9), (4, 7, 8, 9), (4, 8, 5, 9)),
    (5, 7): ((5, 4, 8, 7), (5, 8, 9, 7), (5, 9, 6, 7), (5, 6, 4, 7)),
    (6, 8): ((6, 5, 9, 8), (6, 9, 7, 8), (6, 7, 4, 8), (6, 4, 5, 8)),
}

# Hexahedra27: edge nodes 8-19, face centers 20 (bottom 0123), 21 (front
# 0154), 22 (right 1265), 23 (back 2376), 24 (left 3047), 25 (top 4567),
# body center 26. Eight sub-hexahedra, one per corner.
_HEX27_SUB_HEXAHEDRA = (
    (0, 8, 20, 11, 12, 21, 26, 24),
    (8, 1, 9, 20, 21, 13, 22, 26),
    (20, 9, 2, 10, 26, 22, 14, 23),
    (11, 20, 10, 3, 24, 26, 23, 15),
    (12, 21, 26, 24, 4, 16, 25, 19),
    (21, 13, 22, 26, 16, 5, 17, 25),
    (26, 22, 14, 23, 25, 17, 6, 18),
    (24, 26, 23, 15, 19, 25, 18, 7),
)

_UNSUBDIVIDABLE_WARNED = set()


def _WarnOnceUnsubdividable(geometry_type):
    if geometry_type not in _UNSUBDIVIDABLE_WARNED:
        _UNSUBDIVIDABLE_WARNED.add(geometry_type)
        Kratos.Logger.PrintWarning(
            "Tessellation",
            f"Geometry type {geometry_type} lacks the interior nodes needed for a "
            "conforming subdivision; reducing to corner nodes instead (mid-side "
            "nodes and their field values are dropped).")


def _TessellateHexahedronCorners(corners, mode):
    if mode == "smallest_id_diagonal":
        return _TessellateHexahedronSmallestId(corners)
    return [tuple(corners[i] for i in simplex) for simplex in _HEXAHEDRON_FAN_TABLE]


def _SubdivideTetrahedron10(node_ids, node_coordinates, mode):
    simplices = [tuple(node_ids[i] for i in tet) for tet in _TET10_CORNER_TETS]
    if mode == "smallest_id_diagonal":
        smallest = min(range(4, 10), key=lambda i: node_ids[i])
        diagonal = next(d for d in _TET10_OCTAHEDRON_DIAGONALS if smallest in d)
    else:  # fan: shortest octahedron diagonal
        diagonal = min(
            _TET10_OCTAHEDRON_DIAGONALS,
            key=lambda d: _SquaredDistance(node_coordinates, node_ids[d[0]], node_ids[d[1]]))
    simplices += [tuple(node_ids[i] for i in tet) for tet in _TET10_OCTAHEDRON_TABLES[diagonal]]
    return simplices


def _SubdivideHigherOrder(geometry_type, node_ids, node_coordinates, mode):
    """Subdivision through real mid-side/interior nodes; None when the type
    has no conforming subdivision (linear types and serendipity fallbacks)."""
    if geometry_type in _TRIANGLE6_TYPES:
        return [tuple(node_ids[i] for i in triangle) for triangle in _TRIANGLE6_SUBDIVISION]
    if geometry_type in _QUAD9_TYPES:
        simplices = []
        for sub_quad in _QUAD9_SUB_QUADS:
            simplices += _SplitQuadrilateral(
                tuple(node_ids[i] for i in sub_quad), node_coordinates, mode)
        return simplices
    if geometry_type in _QUAD8_TYPES:
        simplices = [tuple(node_ids[i] for i in triangle) for triangle in _QUAD8_CORNER_TRIANGLES]
        simplices += _SplitQuadrilateral(
            tuple(node_ids[i] for i in _QUAD8_INNER_QUAD), node_coordinates, mode)
        return simplices
    if geometry_type == _GEOMETRY_TYPE.Kratos_Tetrahedra3D10:
        return _SubdivideTetrahedron10(node_ids, node_coordinates, mode)
    if geometry_type == _GEOMETRY_TYPE.Kratos_Hexahedra3D27:
        simplices = []
        for sub_hexahedron in _HEX27_SUB_HEXAHEDRA:
            simplices += _TessellateHexahedronCorners(
                tuple(node_ids[i] for i in sub_hexahedron), mode)
        return simplices
    if geometry_type in _UNSUBDIVIDABLE_TYPES:
        _WarnOnceUnsubdividable(geometry_type)
    return None


# --- Public API -------------------------------------------------------------


def GetSupportedGeometryTypes():
    return frozenset(_CORNER_COUNT.keys())


def TessellateEntity(geometry_type, corner_node_ids, node_coordinates,
                     mode="smallest_id_diagonal", higher_order_mode="reduce"):
    """Decomposes one entity into simplices of Kratos node ids.

    Args:
        geometry_type: The entity's Kratos.GeometryData.KratosGeometryType.
        corner_node_ids: The entity's node ids (higher-order nodes included;
            only the leading corner nodes are used unless
            higher_order_mode == "subdivide").
        node_coordinates: dict-like mapping node id -> numpy (3,) coordinates,
            used for the shortest-diagonal choices of "fan" mode.
        mode: "smallest_id_diagonal" (default, globally consistent faces) or
            "fan" (legacy fixed tables). See the module docstring.
        higher_order_mode: "reduce" (default) or "subdivide". See the module
            docstring.

    Returns:
        list of tuples of Kratos node ids, each a simplex (3 ids for surface
        geometries, 4 for volume geometries).
    """
    if mode not in _TESSELLATION_MODES:
        raise ValueError(f"Unsupported tessellation mode \"{mode}\". Use one of {_TESSELLATION_MODES}.")
    if higher_order_mode not in _HIGHER_ORDER_MODES:
        raise ValueError(
            f"Unsupported higher-order mode \"{higher_order_mode}\". Use one of {_HIGHER_ORDER_MODES}.")
    if geometry_type not in _CORNER_COUNT:
        raise RuntimeError(
            f"Unsupported geometry type for tessellation: {geometry_type}. "
            f"Supported types: {sorted(str(t) for t in _CORNER_COUNT)}.")

    node_ids = tuple(corner_node_ids)

    if higher_order_mode == "subdivide":
        simplices = _SubdivideHigherOrder(geometry_type, node_ids, node_coordinates, mode)
        if simplices is not None:
            return simplices

    corners = node_ids[:_CORNER_COUNT[geometry_type]]

    if geometry_type in _TRIANGLE_TYPES or geometry_type in _TETRAHEDRON_TYPES:
        return [corners]

    if geometry_type in _QUADRILATERAL_TYPES:
        return _SplitQuadrilateral(corners, node_coordinates, mode)

    if geometry_type in _HEXAHEDRON_TYPES:
        return _TessellateHexahedronCorners(corners, mode)

    if geometry_type in _PRISM_TYPES:
        if mode == "smallest_id_diagonal":
            return _TessellatePrismSmallestId(corners)
        return [tuple(corners[i] for i in simplex) for simplex in _PRISM_FAN_TABLE]

    if geometry_type in _PYRAMID_TYPES:
        if mode == "smallest_id_diagonal":
            return _TessellatePyramidSmallestId(corners)
        return [tuple(corners[i] for i in simplex) for simplex in _PYRAMID_FAN_TABLE]

    raise RuntimeError(f"Unhandled geometry type: {geometry_type}.")  # pragma: no cover


def TessellateContainer(container, node_coordinates,
                        mode="smallest_id_diagonal", higher_order_mode="reduce"):
    """Tessellates every entity of a container (Elements or Conditions).

    Args:
        container: Kratos Elements or Conditions container.
        node_coordinates: dict mapping node id -> numpy (3,) coordinates.
        mode: tessellation mode, see TessellateEntity.
        higher_order_mode: higher-order handling, see TessellateEntity.

    Returns:
        (simplex_node_ids, cell_provenance): simplex_node_ids is a list of
        node-id tuples; cell_provenance is a numpy int64 array of shape
        (n_simplices, 2) holding (source_entity_id, sub_cell_index).
    """
    simplex_node_ids = []
    provenance = []
    for entity in container:
        geometry = entity.GetGeometry()
        node_ids = [node.Id for node in geometry]
        simplices = TessellateEntity(
            geometry.GetGeometryType(), node_ids, node_coordinates, mode, higher_order_mode)
        for sub_cell_index, simplex in enumerate(simplices):
            simplex_node_ids.append(simplex)
            provenance.append((entity.Id, sub_cell_index))

    cell_provenance = numpy.array(provenance, dtype=numpy.int64).reshape(-1, 2)
    return simplex_node_ids, cell_provenance
