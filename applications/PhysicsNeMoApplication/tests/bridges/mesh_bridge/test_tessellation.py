import itertools

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import tessellation
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

_GT = Kratos.GeometryData.KratosGeometryType

_UNIT_HEXAHEDRON = {
    0: numpy.array([0.0, 0.0, 0.0]), 1: numpy.array([1.0, 0.0, 0.0]),
    2: numpy.array([1.0, 1.0, 0.0]), 3: numpy.array([0.0, 1.0, 0.0]),
    4: numpy.array([0.0, 0.0, 1.0]), 5: numpy.array([1.0, 0.0, 1.0]),
    6: numpy.array([1.0, 1.0, 1.0]), 7: numpy.array([0.0, 1.0, 1.0]),
}


def _TetrahedronVolume(coordinates, simplex):
    a, b, c, d = (coordinates[node_id] for node_id in simplex)
    return numpy.linalg.det(numpy.array([b - a, c - a, d - a])) / 6.0


def _TriangleArea(coordinates, simplex):
    a, b, c = (coordinates[node_id] for node_id in simplex)
    return 0.5 * numpy.linalg.norm(numpy.cross(b - a, c - a))


def _CreateHexahedraModelPart(model, nx=2):
    """nx unit hexahedra stacked along x."""
    model_part = model.CreateModelPart("hexes")
    props = model_part.CreateNewProperties(1)
    node_id = 0
    for i in range(nx + 1):
        for (y, z) in ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)):
            node_id += 1
            model_part.CreateNewNode(node_id, float(i), y, z)
    for e in range(nx):
        base = 4 * e
        # bottom face (z=0 plane nodes 1,2 y-edge...) -- Kratos hex ordering:
        # 4 nodes of face x=e, then 4 nodes of face x=e+1
        model_part.CreateNewElement(
            "Element3D8N", e + 1,
            [base + 1, base + 5, base + 6, base + 2, base + 4, base + 8, base + 7, base + 3],
            props)
    return model_part


class TestTessellationTables(KratosUnittest.TestCase):
    def test_TriangleIdentity(self):
        simplices = tessellation.TessellateEntity(_GT.Kratos_Triangle2D3, [7, 8, 9], {})
        self.assertEqual(simplices, [(7, 8, 9)])

    def test_TetrahedronIdentity(self):
        simplices = tessellation.TessellateEntity(_GT.Kratos_Tetrahedra3D4, [1, 2, 3, 4], {})
        self.assertEqual(simplices, [(1, 2, 3, 4)])

    def test_QuadrilateralShortestDiagonalFanMode(self):
        # Stretched quad: diagonal 1-3 (indices 0-2) is longer than 2-4 (1-3).
        coordinates = {
            1: numpy.array([0.0, 0.0, 0.0]),
            2: numpy.array([2.0, 0.0, 0.0]),
            3: numpy.array([2.0, 1.0, 0.0]),
            4: numpy.array([0.0, 1.0, 0.0]),
        }
        simplices = tessellation.TessellateEntity(
            _GT.Kratos_Quadrilateral2D4, [1, 2, 3, 4], coordinates, mode="fan")
        self.assertEqual(len(simplices), 2)
        # diag(1,3) has length sqrt(5) == diag(2,4) -> tie goes to 0-2 diagonal
        self.assertIn((1, 2, 3), simplices)
        self.assertIn((1, 3, 4), simplices)

    def test_QuadrilateralSmallestIdDiagonal(self):
        # Smallest id 3 sits at local corner 1 -> diagonal through corners 1-3.
        simplices = tessellation.TessellateEntity(
            _GT.Kratos_Quadrilateral2D4, [9, 3, 7, 5], {})
        self.assertEqual(len(simplices), 2)
        self.assertIn((9, 3, 5), simplices)
        self.assertIn((3, 7, 5), simplices)

    def test_HexahedronFanMode(self):
        simplices = tessellation.TessellateEntity(
            _GT.Kratos_Hexahedra3D8, [1, 2, 3, 4, 5, 6, 7, 8], {}, mode="fan")
        self.assertEqual(len(simplices), 6)
        for simplex in simplices:
            self.assertEqual(len(simplex), 4)
            # fan around the 1-7 diagonal: both endpoints in every tet
            self.assertIn(1, simplex)
            self.assertIn(7, simplex)

    def test_HexahedronSmallestIdMode(self):
        simplices = tessellation.TessellateEntity(
            _GT.Kratos_Hexahedra3D8, [1, 2, 3, 4, 5, 6, 7, 8], {})
        self.assertIn(len(simplices), (5, 6))
        coordinates = {node_id + 1: xyz for node_id, xyz in _UNIT_HEXAHEDRON.items()}
        volumes = [_TetrahedronVolume(coordinates, s) for s in simplices]
        self.assertTrue(all(v > 0.0 for v in volumes))
        self.assertAlmostEqual(sum(volumes), 1.0, places=12)

    def test_PrismDecomposition(self):
        for mode in ("smallest_id_diagonal", "fan"):
            simplices = tessellation.TessellateEntity(
                _GT.Kratos_Prism3D6, [1, 2, 3, 4, 5, 6], {}, mode=mode)
            self.assertEqual(len(simplices), 3)

    def test_PyramidDecomposition(self):
        for mode in ("smallest_id_diagonal", "fan"):
            simplices = tessellation.TessellateEntity(
                _GT.Kratos_Pyramid3D5, [1, 2, 3, 4, 5], {}, mode=mode)
            self.assertEqual(len(simplices), 2)
            for simplex in simplices:
                self.assertIn(5, simplex)  # apex in every tet

    def test_HigherOrderReducesToCorners(self):
        # Tet10: only the 4 corner nodes survive in the default "reduce" mode.
        simplices = tessellation.TessellateEntity(
            _GT.Kratos_Tetrahedra3D10, [1, 2, 3, 4, 11, 12, 13, 14, 15, 16], {})
        self.assertEqual(simplices, [(1, 2, 3, 4)])

    def test_UnsupportedGeometryRaises(self):
        with self.assertRaisesRegex(RuntimeError, "Unsupported geometry type"):
            tessellation.TessellateEntity(_GT.Kratos_Sphere3D1, [1], {})

    def test_UnsupportedModeRaises(self):
        with self.assertRaisesRegex(ValueError, "Unsupported tessellation mode"):
            tessellation.TessellateEntity(_GT.Kratos_Triangle2D3, [1, 2, 3], {}, mode="voronoi")
        with self.assertRaisesRegex(ValueError, "Unsupported higher-order mode"):
            tessellation.TessellateEntity(
                _GT.Kratos_Triangle2D3, [1, 2, 3], {}, higher_order_mode="refine")


class TestDompierreTessellation(KratosUnittest.TestCase):
    """Smallest-id diagonal rule: watertight for arbitrary global numberings."""

    def _CheckHexahedron(self, node_ids):
        coordinates = {node_ids[local]: xyz for local, xyz in _UNIT_HEXAHEDRON.items()}
        simplices = tessellation.TessellateEntity(_GT.Kratos_Hexahedra3D8, node_ids, coordinates)
        volumes = [_TetrahedronVolume(coordinates, s) for s in simplices]
        self.assertTrue(all(v > 1e-12 for v in volumes), f"non-positive tet for ids {node_ids}")
        self.assertAlmostEqual(sum(volumes), 1.0, places=12)
        return simplices

    def test_AllDiagonalConfigurations(self):
        # The diagonal count through the corner opposite the smallest node id
        # takes all values 0..3 over id permutations; every one must yield a
        # positive, volume-conserving decomposition of 5 (n=0) or 6 tets.
        observed_counts = set()
        for permutation in itertools.permutations(range(8)):
            node_ids = [permutation[local] + 1 for local in range(8)]
            simplices = tessellation.TessellateEntity(_GT.Kratos_Hexahedra3D8, node_ids, {})
            observed_counts.add(len(simplices))
            if observed_counts == {5, 6}:
                break
        self.assertEqual(observed_counts, {5, 6})

    def test_RandomIdHexahedra(self):
        rng = numpy.random.default_rng(7)
        for _ in range(100):
            node_ids = [int(i) for i in rng.permutation(numpy.arange(1, 1000))[:8]]
            self._CheckHexahedron(node_ids)

    def test_RandomIdBlockIsWatertight(self):
        # 2x2x2 block of unit hexahedra with random global node ids: every
        # interior triangular face must be shared by exactly two tetrahedra,
        # every boundary face by exactly one. The legacy fan does not have
        # this property for non-translational numberings; the smallest-id
        # rule guarantees it.
        rng = numpy.random.default_rng(3)
        n = 2
        for _ in range(20):
            grid_ids = rng.permutation(numpy.arange(1, 500))[:(n + 1) ** 3].reshape(n + 1, n + 1, n + 1)
            coordinates = {}
            for i, j, k in itertools.product(range(n + 1), repeat=3):
                coordinates[int(grid_ids[i, j, k])] = numpy.array([float(i), float(j), float(k)])
            face_count = {}
            total_volume = 0.0
            for i, j, k in itertools.product(range(n), repeat=3):
                corners = [int(grid_ids[i + di, j + dj, k + dk]) for (di, dj, dk) in
                           [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                            (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]]
                for simplex in tessellation.TessellateEntity(_GT.Kratos_Hexahedra3D8, corners, coordinates):
                    volume = _TetrahedronVolume(coordinates, simplex)
                    self.assertGreater(volume, 1e-12)
                    total_volume += volume
                    for face in itertools.combinations(sorted(simplex), 3):
                        face_count[face] = face_count.get(face, 0) + 1
            self.assertAlmostEqual(total_volume, float(n ** 3), places=10)
            for face, count in face_count.items():
                points = numpy.array([coordinates[node_id] for node_id in face])
                on_boundary = any(
                    numpy.allclose(points[:, axis], value)
                    for axis in range(3) for value in (0.0, float(n)))
                self.assertEqual(count, 1 if on_boundary else 2,
                                 f"face {face} counted {count} times")

    def test_QuadConditionMatchesHexahedronFace(self):
        # A quadrilateral surface condition on a hexahedron face must be
        # triangulated identically to the face of the volume tessellation.
        rng = numpy.random.default_rng(11)
        hex_faces_local = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
                           (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
        for _ in range(50):
            node_ids = [int(i) for i in rng.permutation(numpy.arange(1, 1000))[:8]]
            coordinates = {node_ids[local]: xyz for local, xyz in _UNIT_HEXAHEDRON.items()}
            simplices = tessellation.TessellateEntity(_GT.Kratos_Hexahedra3D8, node_ids, coordinates)
            volume_faces = {face for s in simplices for face in itertools.combinations(sorted(s), 3)}
            for face_local in hex_faces_local:
                quad_ids = [node_ids[local] for local in face_local]
                for triangle in tessellation.TessellateEntity(
                        _GT.Kratos_Quadrilateral3D4, quad_ids, coordinates):
                    self.assertIn(tuple(sorted(triangle)), volume_faces)

    def test_PrismAllSmallestIdPositions(self):
        coordinates_local = [
            numpy.array([0.0, 0.0, 0.0]), numpy.array([1.0, 0.0, 0.0]),
            numpy.array([0.0, 1.0, 0.0]), numpy.array([0.0, 0.0, 1.0]),
            numpy.array([1.0, 0.0, 1.0]), numpy.array([0.0, 1.0, 1.0])]
        rng = numpy.random.default_rng(5)
        for smallest_local in range(6):
            for _ in range(20):
                other_ids = [int(i) for i in rng.permutation(numpy.arange(10, 100))[:5]]
                node_ids = other_ids[:smallest_local] + [2] + other_ids[smallest_local:]
                coordinates = {node_ids[local]: coordinates_local[local] for local in range(6)}
                simplices = tessellation.TessellateEntity(_GT.Kratos_Prism3D6, node_ids, coordinates)
                self.assertEqual(len(simplices), 3)
                volumes = [_TetrahedronVolume(coordinates, s) for s in simplices]
                self.assertTrue(all(v > 1e-12 for v in volumes), f"ids {node_ids}: {volumes}")
                self.assertAlmostEqual(sum(volumes), 0.5, places=12)

    def test_PyramidDiagonalRule(self):
        # Smallest base id at local 1 -> base diagonal 1-3.
        simplices = tessellation.TessellateEntity(_GT.Kratos_Pyramid3D5, [9, 2, 7, 5, 30], {})
        self.assertEqual(simplices, [(9, 2, 5, 30), (2, 7, 5, 30)])
        # Smallest base id at local 2 -> base diagonal 0-2.
        simplices = tessellation.TessellateEntity(_GT.Kratos_Pyramid3D5, [9, 4, 2, 5, 30], {})
        self.assertEqual(simplices, [(9, 4, 2, 30), (9, 2, 5, 30)])


class TestHigherOrderSubdivision(KratosUnittest.TestCase):
    """higher_order_mode="subdivide": subdivision through real mid-side nodes."""

    def setUp(self):
        self.model = Kratos.Model()

    def _CheckProvenanceMeasure(self, model_part, expected_cells, measure="Volume"):
        provenance = domain_mesh_builder.BuildProvenance(model_part, higher_order_mode="subdivide")
        self.assertEqual(provenance.number_of_cells, expected_cells)
        kratos_measure = sum(
            getattr(element.GetGeometry(), measure)() for element in model_part.Elements)
        self.assertAlmostEqual(provenance.ComputeSimplexMeasures().sum(), kratos_measure, places=12)
        return provenance

    def test_Triangle6Subdivision(self):
        model_part = self.model.CreateModelPart("tri6")
        props = model_part.CreateNewProperties(1)
        corners = [(0.0, 0.0), (2.0, 0.0), (0.0, 1.0)]
        mids = [(1.0, 0.0), (1.0, 0.5), (0.0, 0.5)]
        for i, (x, y) in enumerate(corners + mids):
            model_part.CreateNewNode(i + 1, x, y, 0.0)
        model_part.CreateNewElement("Element2D6N", 1, [1, 2, 3, 4, 5, 6], props)
        self._CheckProvenanceMeasure(model_part, 4, "Area")

    def test_Quadrilateral9Subdivision(self):
        model_part = self.model.CreateModelPart("quad9")
        props = model_part.CreateNewProperties(1)
        positions = [(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0),        # corners
                     (1.0, 0.0), (2.0, 0.5), (1.0, 1.0), (0.0, 0.5),        # mid-edges
                     (1.0, 0.5)]                                            # center
        for i, (x, y) in enumerate(positions):
            model_part.CreateNewNode(i + 1, x, y, 0.0)
        model_part.CreateNewElement("Element2D9N", 1, list(range(1, 10)), props)
        self._CheckProvenanceMeasure(model_part, 8, "Area")

    def test_Quadrilateral8Subdivision(self):
        model_part = self.model.CreateModelPart("quad8")
        props = model_part.CreateNewProperties(1)
        positions = [(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0),
                     (1.0, 0.0), (2.0, 0.5), (1.0, 1.0), (0.0, 0.5)]
        for i, (x, y) in enumerate(positions):
            model_part.CreateNewNode(i + 1, x, y, 0.0)
        model_part.CreateNewElement("Element2D8N", 1, list(range(1, 9)), props)
        self._CheckProvenanceMeasure(model_part, 6, "Area")

    def test_Tetrahedron10Subdivision(self):
        model_part = self.model.CreateModelPart("tet10")
        props = model_part.CreateNewProperties(1)
        corners = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        edges = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
        positions = list(corners) + [(corners[a] + corners[b]) / 2.0 for a, b in edges]
        for i, xyz in enumerate(positions):
            model_part.CreateNewNode(i + 1, *[float(c) for c in xyz])
        model_part.CreateNewElement("Element3D10N", 1, list(range(1, 11)), props)
        self._CheckProvenanceMeasure(model_part, 8)

    def test_Hexahedron27Subdivision(self):
        model_part = self.model.CreateModelPart("hex27")
        props = model_part.CreateNewProperties(1)
        corners = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
                               [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]])
        edges = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 5), (2, 6), (3, 7),
                 (4, 5), (5, 6), (6, 7), (7, 4)]
        faces = [(0, 1, 2, 3), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7), (4, 5, 6, 7)]
        positions = list(corners)
        positions += [(corners[a] + corners[b]) / 2.0 for a, b in edges]
        positions += [corners[list(face)].mean(axis=0) for face in faces]
        positions += [corners.mean(axis=0)]
        for i, xyz in enumerate(positions):
            model_part.CreateNewNode(i + 1, *[float(c) for c in xyz])
        model_part.CreateNewElement("Element3D27N", 1, list(range(1, 28)), props)
        # 8 sub-hexahedra of 5 or 6 tets each
        provenance = domain_mesh_builder.BuildProvenance(model_part, higher_order_mode="subdivide")
        self.assertGreaterEqual(provenance.number_of_cells, 40)
        self.assertLessEqual(provenance.number_of_cells, 48)
        self.assertAlmostEqual(
            provenance.ComputeSimplexMeasures().sum(),
            model_part.GetElement(1).GetGeometry().Volume(), places=12)
        # every one of the 27 real nodes is a simplex point
        self.assertEqual(provenance.number_of_points, 27)

    def test_Hexahedron20FallsBackToCorners(self):
        model_part = self.model.CreateModelPart("hex20")
        props = model_part.CreateNewProperties(1)
        corners = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
                               [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]])
        edges = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 5), (2, 6), (3, 7),
                 (4, 5), (5, 6), (6, 7), (7, 4)]
        positions = list(corners) + [(corners[a] + corners[b]) / 2.0 for a, b in edges]
        for i, xyz in enumerate(positions):
            model_part.CreateNewNode(i + 1, *[float(c) for c in xyz])
        model_part.CreateNewElement("Element3D20N", 1, list(range(1, 21)), props)
        provenance = domain_mesh_builder.BuildProvenance(model_part, higher_order_mode="subdivide")
        # corner reduction: only the 8 corner nodes are simplex points
        self.assertEqual(provenance.number_of_points, 8)
        self.assertIn(provenance.number_of_cells, (5, 6))
        self.assertAlmostEqual(provenance.ComputeSimplexMeasures().sum(), 1.0, places=12)

    def test_SubdivisionFieldRoundTrip(self):
        # Nodal gather/scatter stays an exact bijection when mid-side nodes
        # become simplex points.
        model_part = self.model.CreateModelPart("tet10_fields")
        model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        props = model_part.CreateNewProperties(1)
        corners = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        edges = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
        positions = list(corners) + [(corners[a] + corners[b]) / 2.0 for a, b in edges]
        for i, xyz in enumerate(positions):
            node = model_part.CreateNewNode(i + 1, *[float(c) for c in xyz])
            node.SetSolutionStepValue(Kratos.TEMPERATURE, float(10 * (i + 1)))
        model_part.CreateNewElement("Element3D10N", 1, list(range(1, 11)), props)

        provenance = domain_mesh_builder.BuildProvenance(model_part, higher_order_mode="subdivide")
        node_ids = [node.Id for node in model_part.Nodes]
        nodal_values = numpy.array(
            [node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in model_part.Nodes])
        gathered = provenance.GatherNodalField(node_ids, nodal_values)
        self.assertEqual(gathered.shape[0], 10)
        scattered = provenance.ScatterNodalField(node_ids, gathered)
        numpy.testing.assert_allclose(scattered, nodal_values)


class TestVectorizedProvenanceFastPath(KratosUnittest.TestCase):
    """The vectorized homogeneous-simplex path must be bit-identical to the
    general per-entity tessellation path."""

    def _AssertSameProvenance(self, model_part, **kwargs):
        fast = domain_mesh_builder.BuildProvenance(model_part, **kwargs)
        original = domain_mesh_builder._TryVectorizedSimplexProvenance
        domain_mesh_builder._TryVectorizedSimplexProvenance = lambda *a, **k: None
        try:
            general = domain_mesh_builder.BuildProvenance(model_part, **kwargs)
        finally:
            domain_mesh_builder._TryVectorizedSimplexProvenance = original
        numpy.testing.assert_array_equal(fast.simplex_cells, general.simplex_cells)
        numpy.testing.assert_array_equal(fast.point_provenance, general.point_provenance)
        numpy.testing.assert_array_equal(fast.cell_provenance, general.cell_provenance)
        numpy.testing.assert_allclose(fast.simplex_points, general.simplex_points)

    def test_LinearTetrahedraMatchGeneralPath(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("tets")
        props = model_part.CreateNewProperties(1)
        # non-contiguous, permutation-tempting node ids
        for node_id, xyz in [(10, (0, 0, 0)), (2, (1, 0, 0)), (7, (0, 1, 0)),
                             (5, (0, 0, 1)), (12, (1, 1, 1))]:
            model_part.CreateNewNode(node_id, *[float(c) for c in xyz])
        model_part.CreateNewElement("Element3D4N", 3, [10, 2, 7, 5], props)
        model_part.CreateNewElement("Element3D4N", 1, [2, 7, 5, 12], props)
        self._AssertSameProvenance(model_part)

    def test_QuadraticTetrahedraReduceModeMatchGeneralPath(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("tet10")
        props = model_part.CreateNewProperties(1)
        corners = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        edges = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
        positions = list(corners) + [(corners[a] + corners[b]) / 2.0 for a, b in edges]
        for i, xyz in enumerate(positions):
            model_part.CreateNewNode(i + 1, *[float(c) for c in xyz])
        model_part.CreateNewElement("Element3D10N", 1, list(range(1, 11)), props)
        self._AssertSameProvenance(model_part)  # reduce mode: fast path applies
        # subdivide mode must NOT take the fast path (mid-side nodes survive)
        provenance = domain_mesh_builder.BuildProvenance(model_part, higher_order_mode="subdivide")
        self.assertEqual(provenance.number_of_points, 10)

    def test_TriangleConditionsMatchGeneralPath(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("skin")
        props = model_part.CreateNewProperties(1)
        for i, xyz in enumerate([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)]):
            model_part.CreateNewNode(i + 1, *xyz)
        model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1, 2, 3], props)
        model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [1, 3, 4], props)
        self._AssertSameProvenance(model_part, source_container="Conditions")

    def test_MixedContainerFallsBackToGeneralPath(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("mixed")
        props = model_part.CreateNewProperties(1)
        for i, xyz in enumerate([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0),
                                 (0.0, 1.0, 0.0), (0.5, 0.5, 1.0)]):
            model_part.CreateNewNode(i + 1, *xyz)
        model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 5], props)
        model_part.CreateNewElement("Element3D5N", 2, [1, 2, 3, 4, 5], props)
        provenance = domain_mesh_builder.BuildProvenance(model_part)
        self.assertEqual(provenance.number_of_cells, 3)  # 1 tet + pyramid -> 2 tets


class TestTessellationVolumeConservation(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()

    def test_HexahedraVolumeConserved(self):
        model_part = _CreateHexahedraModelPart(self.model, nx=2)
        kratos_volume = sum(element.GetGeometry().Volume() for element in model_part.Elements)

        provenance = domain_mesh_builder.BuildProvenance(model_part)
        self.assertGreaterEqual(provenance.number_of_cells, 10)  # 2 hexes x 5-6 tets
        self.assertLessEqual(provenance.number_of_cells, 12)
        self.assertAlmostEqual(provenance.ComputeSimplexMeasures().sum(), kratos_volume, places=12)

        fan_provenance = domain_mesh_builder.BuildProvenance(model_part, tessellation_mode="fan")
        self.assertEqual(fan_provenance.number_of_cells, 12)  # 2 hexes x 6 tets
        self.assertAlmostEqual(fan_provenance.ComputeSimplexMeasures().sum(), kratos_volume, places=12)

    def test_PrismVolumeConserved(self):
        model_part = self.model.CreateModelPart("prisms")
        props = model_part.CreateNewProperties(1)
        # Two stacked wedges: bottom triangle at z=0/1/2 planes.
        for layer in range(3):
            base = 3 * layer
            model_part.CreateNewNode(base + 1, 0.0, 0.0, float(layer))
            model_part.CreateNewNode(base + 2, 1.0, 0.0, float(layer))
            model_part.CreateNewNode(base + 3, 0.0, 1.0, float(layer))
        model_part.CreateNewElement("Element3D6N", 1, [1, 2, 3, 4, 5, 6], props)
        model_part.CreateNewElement("Element3D6N", 2, [4, 5, 6, 7, 8, 9], props)

        provenance = domain_mesh_builder.BuildProvenance(model_part)
        self.assertEqual(provenance.number_of_cells, 6)  # 2 prisms x 3 tets
        kratos_volume = sum(element.GetGeometry().Volume() for element in model_part.Elements)
        self.assertAlmostEqual(provenance.ComputeSimplexMeasures().sum(), kratos_volume, places=12)

    def test_PyramidVolumeConserved(self):
        model_part = self.model.CreateModelPart("pyramids")
        props = model_part.CreateNewProperties(1)
        for i, xyz in enumerate([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0),
                                 (0.0, 1.0, 0.0), (0.5, 0.5, 1.0)]):
            model_part.CreateNewNode(i + 1, *xyz)
        model_part.CreateNewElement("Element3D5N", 1, [1, 2, 3, 4, 5], props)

        provenance = domain_mesh_builder.BuildProvenance(model_part)
        self.assertEqual(provenance.number_of_cells, 2)
        kratos_volume = model_part.GetElement(1).GetGeometry().Volume()
        self.assertAlmostEqual(provenance.ComputeSimplexMeasures().sum(), kratos_volume, places=12)

    def test_SurfaceConditionsAreaConserved(self):
        model_part = self.model.CreateModelPart("surface")
        props = model_part.CreateNewProperties(1)
        # Two unit quads in the z=0 plane, tessellated from Conditions.
        for i, xyz in enumerate([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                                 (0.0, 1.0, 0.0), (1.0, 1.0, 0.0), (2.0, 1.0, 0.0)]):
            model_part.CreateNewNode(i + 1, *xyz)
        model_part.CreateNewCondition("SurfaceCondition3D4N", 1, [1, 2, 5, 4], props)
        model_part.CreateNewCondition("SurfaceCondition3D4N", 2, [2, 3, 6, 5], props)

        provenance = domain_mesh_builder.BuildProvenance(model_part, "Conditions")
        self.assertEqual(provenance.number_of_cells, 4)  # 2 quads x 2 triangles
        self.assertEqual(provenance.simplex_cells.shape[1], 3)
        kratos_area = sum(condition.GetGeometry().Area() for condition in model_part.Conditions)
        self.assertAlmostEqual(provenance.ComputeSimplexMeasures().sum(), kratos_area, places=12)

    def test_SharedFaceConsistency(self):
        # The interface face between two neighbouring hexahedra must be
        # triangulated identically by both, in every mode (this numbering is
        # translationally consistent, so even the fan complies).
        for mode in ("smallest_id_diagonal", "fan"):
            model = Kratos.Model()
            model_part = _CreateHexahedraModelPart(model, nx=2)
            provenance = domain_mesh_builder.BuildProvenance(model_part, tessellation_mode=mode)

            interface_node_ids = {5, 6, 7, 8}  # nodes at x=1
            interface_point_indices = set(provenance.GetPointIndices(sorted(interface_node_ids)).tolist())

            def interface_triangles(element_row):
                mask = provenance.cell_provenance[:, 0] == element_row
                triangles = set()
                for cell in provenance.simplex_cells[mask]:
                    face = frozenset(int(p) for p in cell if int(p) in interface_point_indices)
                    if len(face) == 3:
                        triangles.add(face)
                return triangles

            self.assertEqual(interface_triangles(1), interface_triangles(2))
            self.assertEqual(len(interface_triangles(1)), 2)


if __name__ == '__main__':
    KratosUnittest.main()
