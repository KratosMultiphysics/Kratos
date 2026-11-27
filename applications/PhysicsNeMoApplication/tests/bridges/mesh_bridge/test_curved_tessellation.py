import itertools

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import curved_tessellation
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

_GT = Kratos.GeometryData.KratosGeometryType

# nodal local coordinates per curved type (Kratos orderings)
_TRI6_LOCALS = [(0, 0), (1, 0), (0, 1), (0.5, 0), (0.5, 0.5), (0, 0.5)]
_QUAD_CORNERS = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
_QUAD_MIDS = [(0, -1), (1, 0), (0, 1), (-1, 0)]
_QUAD8_LOCALS = _QUAD_CORNERS + _QUAD_MIDS
_QUAD9_LOCALS = _QUAD8_LOCALS + [(0, 0)]
_TET10_LOCALS = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
                 (0.5, 0, 0), (0.5, 0.5, 0), (0, 0.5, 0),
                 (0, 0, 0.5), (0.5, 0, 0.5), (0, 0.5, 0.5)]
_FACTOR_LOCAL = {0: -1.0, 1: 1.0, 2: 0.0}
_HEX27_LOCALS = [tuple(_FACTOR_LOCAL[i] for i in ijk)
                 for ijk in curved_tessellation._HEX27_TENSOR_ORDER]

_ELEMENT_NAMES = {
    _GT.Kratos_Triangle2D6: ("Element2D6N", _TRI6_LOCALS),
    _GT.Kratos_Quadrilateral2D8: ("Element2D8N", _QUAD8_LOCALS),
    _GT.Kratos_Quadrilateral2D9: ("Element2D9N", _QUAD9_LOCALS),
    _GT.Kratos_Tetrahedra3D10: ("Element3D10N", _TET10_LOCALS),
    _GT.Kratos_Hexahedra3D27: ("Element3D27N", _HEX27_LOCALS),
}


def _DistortedPositions(locals_list, seed=3, amplitude=0.08):
    """Node positions = local coords (padded to 3D) + a random distortion."""
    rng = numpy.random.default_rng(seed)
    positions = []
    for local in locals_list:
        point = numpy.zeros(3)
        point[:len(local)] = local
        positions.append(point + rng.uniform(-amplitude, amplitude, 3))
    return positions


def _CreateSingleElementPart(model, geometry_type, positions, name="Main"):
    model_part = model.CreateModelPart(name)
    model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
    model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
    properties = model_part.CreateNewProperties(1)
    for i, xyz in enumerate(positions):
        model_part.CreateNewNode(i + 1, *[float(c) for c in xyz])
    element_name, _ = _ELEMENT_NAMES[geometry_type]
    model_part.CreateNewElement(element_name, 1, list(range(1, len(positions) + 1)), properties)
    return model_part


def _Hex27Positions(x_range, distortion=None):
    """27 node positions of a hex spanning x_range x [0,1] x [0,1]."""
    positions = []
    for local in _HEX27_LOCALS:
        point = numpy.array([
            x_range[0] + (local[0] + 1.0) / 2.0 * (x_range[1] - x_range[0]),
            (local[1] + 1.0) / 2.0, (local[2] + 1.0) / 2.0])
        if distortion is not None:
            point = distortion(point)
        positions.append(point)
    return positions


def _FaceCounts(provenance):
    counts = {}
    for simplex in provenance.simplex_cells:
        for face in itertools.combinations(sorted(int(v) for v in simplex), 3):
            counts[face] = counts.get(face, 0) + 1
    return counts


def _SimplexIdSets(provenance):
    """Simplices as frozensets of node ids (only valid without synthetics)."""
    ids = provenance.point_provenance
    return {frozenset(int(ids[v]) for v in simplex) for simplex in provenance.simplex_cells}


class TestCurvedShapeFunctions(KratosUnittest.TestCase):
    def _CheckType(self, geometry_type, locals_list, places=12):
        rng = numpy.random.default_rng(11)
        dimension = len(locals_list[0])
        if geometry_type in (_GT.Kratos_Triangle2D6,):
            a = rng.uniform(0.0, 1.0, (20, 1))
            b = rng.uniform(0.0, 1.0, (20, 1)) * (1.0 - a)
            samples = numpy.hstack([a, b])
        elif geometry_type == _GT.Kratos_Tetrahedra3D10:
            raw = rng.uniform(0.0, 1.0, (20, 3))
            samples = raw / numpy.maximum(raw.sum(axis=1), 1.0)[:, None] * 0.9
        else:
            samples = rng.uniform(-1.0, 1.0, (20, dimension))

        weights = curved_tessellation.EvaluateShapeFunctions(geometry_type, samples)
        # partition of unity
        numpy.testing.assert_allclose(weights.sum(axis=1), 1.0, atol=1e-12)
        # nodal delta property
        nodal = curved_tessellation.EvaluateShapeFunctions(geometry_type, numpy.array(locals_list))
        numpy.testing.assert_allclose(nodal, numpy.eye(len(locals_list)), atol=1e-12)
        # exact reproduction of a complete quadratic polynomial
        def quadratic(points):
            padded = numpy.zeros((len(points), 3))
            padded[:, :points.shape[1]] = points
            x, y, z = padded[:, 0], padded[:, 1], padded[:, 2]
            return 1.0 + 2.0 * x - y + 0.5 * z + x * x + 3.0 * y * y + z * z + x * y - y * z
        numpy.testing.assert_allclose(
            weights @ quadratic(numpy.array(locals_list, dtype=float)),
            quadratic(samples), atol=1e-11)
        # cross-check against the Kratos isoparametric map
        model = Kratos.Model()
        model_part = _CreateSingleElementPart(
            model, geometry_type, _DistortedPositions(locals_list))
        geometry = model_part.GetElement(1).GetGeometry()
        coordinates = numpy.array([[node.X, node.Y, node.Z] for node in geometry])
        for local in samples[:5]:
            padded = numpy.zeros(3)
            padded[:len(local)] = local
            kratos_point = numpy.array(geometry.GlobalCoordinates(Kratos.Vector(padded)))
            numpy_point = curved_tessellation.EvaluateShapeFunctions(
                geometry_type, local[None])[0] @ coordinates
            numpy.testing.assert_allclose(numpy_point, kratos_point, atol=1e-12)

    def test_AllCurvedTypes(self):
        for geometry_type, (_, locals_list) in _ELEMENT_NAMES.items():
            with self.subTest(geometry_type=str(geometry_type)):
                self._CheckType(geometry_type, locals_list)


class TestCurvedLevelOneEquivalence(KratosUnittest.TestCase):
    def _Compare(self, geometry_type):
        model = Kratos.Model()
        _, locals_list = _ELEMENT_NAMES[geometry_type]
        model_part = _CreateSingleElementPart(
            model, geometry_type, _DistortedPositions(locals_list))

        curved = domain_mesh_builder.BuildProvenance(
            model_part, higher_order_mode="curved", curved_refinement_levels=1)
        subdivided = domain_mesh_builder.BuildProvenance(
            model_part, higher_order_mode="subdivide")

        self.assertEqual(curved.number_of_synthetic_points, 0)
        numpy.testing.assert_array_equal(curved.point_provenance, subdivided.point_provenance)
        self.assertEqual(_SimplexIdSets(curved), _SimplexIdSets(subdivided))
        self.assertAlmostEqual(
            float(curved.ComputeSimplexMeasures().sum()),
            float(subdivided.ComputeSimplexMeasures().sum()), places=12)

    def test_Quadratic9Equivalence(self):
        self._Compare(_GT.Kratos_Quadrilateral2D9)

    def test_Triangle6Equivalence(self):
        self._Compare(_GT.Kratos_Triangle2D6)

    def test_Tetrahedra10Equivalence(self):
        # sequential node ids make subdivide's smallest-mid-id octahedron rule
        # pick the (4, 9) diagonal - the curved mode's fixed parametric choice
        self._Compare(_GT.Kratos_Tetrahedra3D10)

    def test_Hexahedra27Equivalence(self):
        self._Compare(_GT.Kratos_Hexahedra3D27)

    def test_Quadrilateral8CenterIsSynthetic(self):
        model = Kratos.Model()
        model_part = _CreateSingleElementPart(
            model, _GT.Kratos_Quadrilateral2D8, _DistortedPositions(_QUAD8_LOCALS))
        curved = domain_mesh_builder.BuildProvenance(
            model_part, higher_order_mode="curved", curved_refinement_levels=1)
        self.assertEqual(curved.number_of_synthetic_points, 1)
        self.assertEqual(curved.number_of_points, 9)
        self.assertEqual(len(curved.simplex_cells), 8)  # 4 sub-quads x 2 triangles
        # the serendipity center weights: corners -1/4, mid-sides 1/2
        numpy.testing.assert_allclose(
            curved.synthetic_weights[0][:8],
            [-0.25] * 4 + [0.5] * 4, atol=1e-12)


class TestCurvedMeasureConvergence(KratosUnittest.TestCase):
    def _Measures(self, model_part, levels_list):
        measures = []
        for levels in levels_list:
            provenance = domain_mesh_builder.BuildProvenance(
                model_part, higher_order_mode="curved", curved_refinement_levels=levels)
            measures.append(float(provenance.ComputeSimplexMeasures().sum()))
        subdivided = domain_mesh_builder.BuildProvenance(
            model_part, higher_order_mode="subdivide")
        return measures, float(subdivided.ComputeSimplexMeasures().sum())

    def test_CircularTriangleAreaConverges(self):
        # corners on the unit circle, mid-side node 3 ON the arc: the exact
        # isoparametric area comes from Kratos' quadratic quadrature (exact
        # for the quadratic Jacobian of a Tri6)
        arc = numpy.sqrt(0.5)
        positions = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 0.0),
                     (arc, arc, 0.0), (0.0, 0.5, 0.0), (0.5, 0.0, 0.0)]
        model = Kratos.Model()
        model_part = _CreateSingleElementPart(model, _GT.Kratos_Triangle2D6, positions)
        exact = model_part.GetElement(1).GetGeometry().Area()

        measures, straight = self._Measures(model_part, (1, 2, 3))
        errors = [abs(m - exact) for m in measures]
        self.assertLess(errors[1], errors[0])
        self.assertLess(errors[2], errors[1])
        self.assertLess(errors[2], errors[0] / 4.0)
        # every curved level beats the straight-edged subdivision
        self.assertLess(errors[0], abs(straight - exact) + 1e-14)
        self.assertLess(errors[2], abs(straight - exact))

    def test_CurvedTetrahedronVolumeCauchyConverges(self):
        positions = _DistortedPositions(_TET10_LOCALS, seed=5, amplitude=0.06)
        model = Kratos.Model()
        model_part = _CreateSingleElementPart(model, _GT.Kratos_Tetrahedra3D10, positions)
        (m1, m2, m3), straight = self._Measures(model_part, (1, 2, 3))
        self.assertLess(abs(m3 - m2), abs(m2 - m1))
        self.assertLess(abs(m2 - m3), abs(straight - m3) + 1e-14)

    def test_CurvedHexahedronVolumeCauchyConverges(self):
        def distortion(point):
            # radial bulge: corners fixed, mid nodes inflated outward - a
            # genuine O(a) volume change (shear-type distortions are
            # volume-preserving and would make every measure identical)
            radial = point - 0.5
            inflation = 0.2 * (0.75 - float(radial @ radial)) / 0.75
            return 0.5 + radial * (1.0 + inflation)
        positions = _Hex27Positions((0.0, 1.0), distortion)
        model = Kratos.Model()
        model_part = _CreateSingleElementPart(model, _GT.Kratos_Hexahedra3D27, positions)
        (m1, m2, m3), straight = self._Measures(model_part, (1, 2, 3))
        self.assertLess(abs(m3 - m2), abs(m2 - m1))


class TestCurvedWatertightness(KratosUnittest.TestCase):
    def _CreateTwoHexPart(self, model, permute_seed=None):
        """Two Hex27 sharing the x = 1 face, with distorted shared-face nodes."""
        def distortion(point):
            if abs(point[0] - 1.0) < 1e-12 and 0.0 < point[1] < 1.0 or \
               abs(point[0] - 1.0) < 1e-12 and 0.0 < point[2] < 1.0:
                bump = 0.08 * numpy.sin(numpy.pi * point[1]) * numpy.sin(numpy.pi * point[2])
                return point + numpy.array([bump, 0.0, 0.0])
            return point

        left = _Hex27Positions((0.0, 1.0), distortion)
        right = _Hex27Positions((1.0, 2.0), distortion)

        position_to_id = {}
        assignments = []  # per hex: list of node ids
        all_positions = []
        for positions in (left, right):
            ids = []
            for point in positions:
                key = tuple(numpy.round(point, 9))
                if key not in position_to_id:
                    position_to_id[key] = len(position_to_id) + 1
                    all_positions.append(point)
                ids.append(position_to_id[key])
            assignments.append(ids)

        if permute_seed is not None:
            rng = numpy.random.default_rng(permute_seed)
            permutation = rng.permutation(numpy.arange(1, len(all_positions) + 1))
            remap = {old: int(permutation[old - 1]) for old in range(1, len(all_positions) + 1)}
            assignments = [[remap[i] for i in ids] for ids in assignments]
            reordered = [None] * len(all_positions)
            for old, new in remap.items():
                reordered[new - 1] = all_positions[old - 1]
            all_positions = reordered

        model_part = model.CreateModelPart("TwoHexes")
        properties = model_part.CreateNewProperties(1)
        for i, point in enumerate(all_positions):
            model_part.CreateNewNode(i + 1, *[float(c) for c in point])
        for e, ids in enumerate(assignments):
            model_part.CreateNewElement("Element3D27N", e + 1, ids, properties)
        return model_part

    def _CheckFaceCounts(self, provenance, expected_boundary):
        counts = _FaceCounts(provenance)
        self.assertLessEqual(max(counts.values()), 2)
        boundary = sum(1 for c in counts.values() if c == 1)
        self.assertEqual(boundary, expected_boundary)
        # no coincident duplicate points
        unique_points = numpy.unique(numpy.round(provenance.simplex_points, 9), axis=0)
        self.assertEqual(len(unique_points), provenance.number_of_points)

    def test_TwoCurvedHexahedraAreWatertight(self):
        for permute_seed in (None, 1, 2):
            model = Kratos.Model()
            model_part = self._CreateTwoHexPart(model, permute_seed)
            for levels in (1, 2):
                provenance = domain_mesh_builder.BuildProvenance(
                    model_part, higher_order_mode="curved", curved_refinement_levels=levels)
                # 10 boundary parent faces, each 2 * 4^k triangles
                self._CheckFaceCounts(provenance, 10 * 2 * 4 ** levels)

    def test_TwoCurvedTetrahedraAreWatertight(self):
        # two Tet10 sharing the (1,0,0)-(0,1,0)-(0,0,1) face, curved mids
        base = {(1, 0, 0): 2, (0, 1, 0): 3, (0, 0, 1): 4}
        positions = {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (0.0, 1.0, 0.0),
                     4: (0.0, 0.0, 1.0), 11: (1.0, 1.0, 1.0)}

        def midpoint(a, b, bump=0.0):
            p = (numpy.array(positions[a]) + numpy.array(positions[b])) / 2.0
            return tuple(p + bump)

        # shared-face mids (2-3, 3-4, 2-4) bumped off the plane
        positions[5] = midpoint(1, 2)
        positions[6] = midpoint(2, 3, 0.05)
        positions[7] = midpoint(1, 3)
        positions[8] = midpoint(1, 4)
        positions[9] = midpoint(2, 4, 0.04)
        positions[10] = midpoint(3, 4, 0.03)
        positions[12] = midpoint(2, 11)
        positions[13] = midpoint(3, 11)
        positions[14] = midpoint(4, 11)

        model = Kratos.Model()
        model_part = model.CreateModelPart("TwoTets")
        properties = model_part.CreateNewProperties(1)
        for node_id, xyz in positions.items():
            model_part.CreateNewNode(node_id, *[float(c) for c in xyz])
        # tet 1: corners 1,2,3,4; mids 5(1-2), 6(2-3), 7(1-3), 8(1-4), 9(2-4), 10(3-4)
        model_part.CreateNewElement("Element3D10N", 1, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], properties)
        # tet 2: corners 2,3,4,11; mids 6(2-3), 10(3-4), 9(2-4)... ordering:
        # mids 4=(0,1)->2-3=6, 5=(1,2)->3-4=10, 6=(0,2)->2-4=9,
        # 7=(0,3)->2-11=12, 8=(1,3)->3-11=13, 9=(2,3)->4-11=14
        model_part.CreateNewElement("Element3D10N", 2, [2, 3, 4, 11, 6, 10, 9, 12, 13, 14], properties)

        for levels in (1, 2):
            provenance = domain_mesh_builder.BuildProvenance(
                model_part, higher_order_mode="curved", curved_refinement_levels=levels)
            self._CheckFaceCounts(provenance, 6 * 4 ** levels)

    def test_ConditionMatchesVolumeBoundary(self):
        # a Quadrilateral3D9 condition on the shared 9 nodes of a hex face
        # must triangulate exactly like the volume boundary there
        model = Kratos.Model()
        model_part = self._CreateTwoHexPart(model)
        # the shared face x = 1 carries the 9 nodes common to both hexes
        left_ids = set(node.Id for node in model_part.GetElement(1).GetGeometry())
        right_ids = set(node.Id for node in model_part.GetElement(2).GetGeometry())
        shared = sorted(left_ids & right_ids)
        self.assertEqual(len(shared), 9)
        coordinates = {node.Id: numpy.array([node.X, node.Y, node.Z])
                       for node in model_part.Nodes}
        # face frame: 4 corners (y,z extremes), then mids, then center
        def rank(node_id):
            _, y, z = coordinates[node_id]
            return (round(y, 6), round(z, 6))
        by_pos = {rank(i): i for i in shared}
        cycle = [by_pos[(0.0, 0.0)], by_pos[(1.0, 0.0)], by_pos[(1.0, 1.0)], by_pos[(0.0, 1.0)]]
        mids = [by_pos[(0.5, 0.0)], by_pos[(1.0, 0.5)], by_pos[(0.5, 1.0)], by_pos[(0.0, 0.5)]]
        center = by_pos[(0.5, 0.5)]
        properties = model_part.GetElement(1).Properties
        model_part.CreateNewCondition(
            "SurfaceCondition3D9N", 1, cycle + mids + [center], properties)

        volume = domain_mesh_builder.BuildProvenance(
            model_part, higher_order_mode="curved", curved_refinement_levels=2)
        surface = domain_mesh_builder.BuildProvenance(
            model_part, "Conditions", higher_order_mode="curved", curved_refinement_levels=2)

        def triangle_set(provenance, simplex_size):
            triangles = set()
            for simplex in provenance.simplex_cells:
                faces = ([tuple(simplex)] if simplex_size == 3
                         else itertools.combinations(simplex, 3))
                for face in faces:
                    points = numpy.round(provenance.simplex_points[list(face)], 8)
                    triangles.add(frozenset(map(tuple, points)))
            return triangles

        surface_triangles = triangle_set(surface, 3)
        volume_triangles = triangle_set(volume, 4)
        self.assertTrue(surface_triangles <= volume_triangles,
                        "condition triangulation diverges from the volume boundary")

    def test_MixedTriangleQuadEdgeMerges(self):
        # Tri6 and Quad9 sharing a (curved) edge: the edge lattice merges
        tri_positions = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (0.0, 1.0),
                         4: (0.5, -0.06), 5: (0.5, 0.5), 6: (0.0, 0.5)}
        quad_positions = {2: (1.0, 0.0), 1: (0.0, 0.0),  # shared edge 1-2 (mid 4)
                          7: (0.0, -1.0), 8: (1.0, -1.0),
                          9: (0.5, -1.0), 10: (0.0, -0.5), 11: (1.0, -0.5),
                          12: (0.5, -0.5)}
        model = Kratos.Model()
        model_part = model.CreateModelPart("Mixed")
        properties = model_part.CreateNewProperties(1)
        seen = {}
        for source in (tri_positions, quad_positions):
            for node_id, xy in source.items():
                if node_id not in seen:
                    seen[node_id] = True
                    model_part.CreateNewNode(node_id, float(xy[0]), float(xy[1]), 0.0)
        model_part.CreateNewElement("Element2D6N", 1, [1, 2, 3, 4, 5, 6], properties)
        # quad corners CCW: 7 (bottom-left), 8, 2, 1; mids 9 (7-8), 11 (8-2),
        # 4 (2-1, the SHARED curved mid), 10 (1-7); center 12
        model_part.CreateNewElement("Element2D9N", 2, [7, 8, 2, 1, 9, 11, 4, 10, 12], properties)

        provenance = domain_mesh_builder.BuildProvenance(
            model_part, higher_order_mode="curved", curved_refinement_levels=2)
        unique_points = numpy.unique(numpy.round(provenance.simplex_points, 9), axis=0)
        self.assertEqual(len(unique_points), provenance.number_of_points)
        # tri lattice 15 + quad lattice 25 - 5 shared edge points
        self.assertEqual(provenance.number_of_points, 35)


class TestCurvedProvenanceAndGather(KratosUnittest.TestCase):
    def _AffineTetPart(self, model):
        corners = numpy.array([[0.0, 0.0, 0.0], [1.2, 0.1, 0.0],
                               [0.1, 1.1, 0.0], [-0.1, 0.2, 1.3]])
        edges = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
        positions = list(corners) + [(corners[a] + corners[b]) / 2.0 for a, b in edges]
        return _CreateSingleElementPart(model, _GT.Kratos_Tetrahedra3D10, positions)

    def test_ProvenanceLayout(self):
        model = Kratos.Model()
        model_part = self._AffineTetPart(model)
        provenance = domain_mesh_builder.BuildProvenance(
            model_part, higher_order_mode="curved", curved_refinement_levels=2)

        synthetic_count = provenance.number_of_synthetic_points
        self.assertGreater(synthetic_count, 0)
        real_count = provenance.number_of_points - synthetic_count
        self.assertTrue((provenance.point_provenance[:real_count] >= 0).all())
        self.assertTrue((provenance.point_provenance[real_count:] == -1).all())
        self.assertTrue((numpy.diff(provenance.point_provenance[:real_count]) > 0).all())
        numpy.testing.assert_allclose(
            provenance.synthetic_weights.sum(axis=1), 1.0, atol=1e-12)
        self.assertTrue((provenance.synthetic_parent_ids == 1).all())
        self.assertNotIn(-1, provenance._node_id_to_point)

    def test_QuadraticFieldGatherIsExactAtSyntheticPoints(self):
        def field(points):
            x, y, z = points[:, 0], points[:, 1], points[:, 2]
            return x * x + 2.0 * y * y + 3.0 * z * z + x * y - z

        for build in (self._AffineTetPart,
                      lambda model: _CreateSingleElementPart(
                          model, _GT.Kratos_Hexahedra3D27, _Hex27Positions((0.0, 1.0)))):
            model = Kratos.Model()
            model_part = build(model)
            provenance = domain_mesh_builder.BuildProvenance(
                model_part, higher_order_mode="curved", curved_refinement_levels=2)

            node_ids = [node.Id for node in model_part.Nodes]
            coordinates = numpy.array([[n.X, n.Y, n.Z] for n in model_part.Nodes])
            nodal = field(coordinates)
            gathered = provenance.GatherNodalField(node_ids, nodal)
            # affine map -> quadratic fields interpolate exactly, synthetic included
            numpy.testing.assert_allclose(
                gathered, field(provenance.simplex_points), atol=1e-10)

            # round trip stays exact on real nodes
            scattered = provenance.ScatterNodalField(node_ids, gathered)
            numpy.testing.assert_allclose(scattered, nodal, atol=1e-12)

    def test_VectorFieldGather(self):
        model = Kratos.Model()
        model_part = self._AffineTetPart(model)
        provenance = domain_mesh_builder.BuildProvenance(
            model_part, higher_order_mode="curved", curved_refinement_levels=2)
        node_ids = [node.Id for node in model_part.Nodes]
        nodal = numpy.array([[n.X, 2.0 * n.Y, n.X + n.Z] for n in model_part.Nodes])
        gathered = provenance.GatherNodalField(node_ids, nodal)
        expected = numpy.stack([provenance.simplex_points[:, 0],
                                2.0 * provenance.simplex_points[:, 1],
                                provenance.simplex_points[:, 0] + provenance.simplex_points[:, 2]],
                               axis=1)
        numpy.testing.assert_allclose(gathered, expected, atol=1e-10)

    def test_ScatterFieldBackWithCurvedProvenance(self):
        model = Kratos.Model()
        model_part = self._AffineTetPart(model)
        provenance = domain_mesh_builder.BuildProvenance(
            model_part, higher_order_mode="curved", curved_refinement_levels=2)
        node_ids = [node.Id for node in model_part.Nodes]
        for node in model_part.Nodes:
            node.SetSolutionStepValue(Kratos.TEMPERATURE, 1.0 + node.X)
        nodal = numpy.array([1.0 + n.X for n in model_part.Nodes])
        prediction = provenance.GatherNodalField(node_ids, nodal)

        domain_mesh_builder.ScatterFieldBack(
            provenance, prediction, model_part, Kratos.TEMPERATURE, "node_historical")
        for node in model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE), 1.0 + node.X, places=10)


class TestCurvedValidation(KratosUnittest.TestCase):
    def test_FastPathDoesNotFireForCurvedQuadratics(self):
        model = Kratos.Model()
        corners = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        edges = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
        positions = list(corners) + [(corners[a] + corners[b]) / 2.0 for a, b in edges]
        model_part = _CreateSingleElementPart(model, _GT.Kratos_Tetrahedra3D10, positions)
        provenance = domain_mesh_builder.BuildProvenance(
            model_part, higher_order_mode="curved", curved_refinement_levels=2)
        self.assertGreater(provenance.number_of_synthetic_points, 0)

    def test_CurvedOnLinearTypesIsANoOp(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("LinearTets")
        properties = model_part.CreateNewProperties(1)
        for i, xyz in enumerate([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)]):
            model_part.CreateNewNode(i + 1, *[float(c) for c in xyz])
        model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], properties)
        model_part.CreateNewElement("Element3D4N", 2, [2, 3, 4, 5], properties)

        curved = domain_mesh_builder.BuildProvenance(model_part, higher_order_mode="curved")
        default = domain_mesh_builder.BuildProvenance(model_part)
        self.assertEqual(curved.number_of_synthetic_points, 0)
        numpy.testing.assert_array_equal(curved.point_provenance, default.point_provenance)
        numpy.testing.assert_array_equal(curved.simplex_cells, default.simplex_cells)
        numpy.testing.assert_array_equal(curved.cell_provenance, default.cell_provenance)

    def test_InvalidCombinationsRaise(self):
        model = Kratos.Model()
        model_part = _CreateSingleElementPart(
            model, _GT.Kratos_Triangle2D6, _DistortedPositions(_TRI6_LOCALS))
        with self.assertRaisesRegex(ValueError, "smallest_id_diagonal"):
            domain_mesh_builder.BuildProvenance(
                model_part, tessellation_mode="fan", higher_order_mode="curved")
        with self.assertRaisesRegex(ValueError, "curved_refinement_levels"):
            domain_mesh_builder.BuildProvenance(
                model_part, higher_order_mode="curved", curved_refinement_levels=0)


if __name__ == '__main__':
    KratosUnittest.main()
