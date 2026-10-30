import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
# Kratos hexahedron ordering -> 6 tetrahedra (fan around the 0-6 diagonal),
# reused from the mesh bridge's decomposition table.
_HEX_TO_TETS = ((0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6))


def CreateStructuredTetModelPart(model, name, divisions, historical_variables=(Kratos.PRESSURE,),
                                 extent=(1.0, 1.0, 1.0)):
    """Structured tetrahedral box mesh with `divisions` cells per axis.

    Defaults to the unit cube, which is what every existing caller wants.
    Pass a non-uniform `extent` to get an ANISOTROPIC box: on [0,1]^3 every
    length scale is 1, so a per-axis or 1/L factor is the identity and a
    length-scale bug is invisible - which is how a rescaled-PDE bug once
    survived a full test suite.
    """
    model_part = model.CreateModelPart(name)
    for variable in historical_variables:
        model_part.AddNodalSolutionStepVariable(variable)
    props = model_part.CreateNewProperties(1)

    n = divisions + 1
    node_id = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                node_id += 1
                model_part.CreateNewNode(node_id,
                                         extent[0] * i / divisions,
                                         extent[1] * j / divisions,
                                         extent[2] * k / divisions)

    def nid(i, j, k):
        return i * n * n + j * n + k + 1

    element_id = 0
    for i in range(divisions):
        for j in range(divisions):
            for k in range(divisions):
                # Kratos hex ordering: bottom face CCW, then top face CCW.
                corners = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                           nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                for tet in _HEX_TO_TETS:
                    element_id += 1
                    model_part.CreateNewElement("Element3D4N", element_id, [corners[c] for c in tet], props)
    return model_part


def _SetLinearField(model_part, variable=Kratos.PRESSURE):
    """u(x,y,z) = 1 + 2x + 3y - z: exactly representable by linear tets."""
    for node in model_part.Nodes:
        node.SetSolutionStepValue(variable, 1.0 + 2.0 * node.X + 3.0 * node.Y - node.Z)


def _ExactLinear(points):
    points = numpy.asarray(points)
    return 1.0 + 2.0 * points[:, 0] + 3.0 * points[:, 1] - points[:, 2]


class TestGridBridge(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(self.model, "Main", divisions=3)
        _SetLinearField(self.model_part)

    def test_BoundingBox(self):
        low, high = grid_bridge.ComputeBoundingBox(self.model_part)
        self.assertTrue(numpy.allclose(low, [0.0, 0.0, 0.0]))
        self.assertTrue(numpy.allclose(high, [1.0, 1.0, 1.0]))

    def test_LinearFieldSampledExactly(self):
        grid, box = grid_bridge.SampleFieldsOnGrid(
            self.model_part, [("PRESSURE", "node_historical")], (5, 4, 6))
        self.assertEqual(grid.shape, (1, 5, 4, 6))
        points = grid_bridge._GridPointCoordinates((5, 4, 6), box)
        self.assertTrue(numpy.allclose(grid.reshape(1, -1)[0], _ExactLinear(points), atol=1e-12))

    def test_OutOfMeshFillValue(self):
        box = (numpy.array([-1.0, -1.0, -1.0]), numpy.array([2.0, 2.0, 2.0]))
        grid, _ = grid_bridge.SampleFieldsOnGrid(
            self.model_part, [("PRESSURE", "node_historical")], (4, 4, 4),
            bounding_box=box, fill_value=-999.0)
        self.assertEqual(grid[0, 0, 0, 0], -999.0)   # corner far outside the cube
        self.assertNotEqual(grid[0, 1, 1, 1], -999.0)  # interior lattice point inside

    def test_TrilinearInterpolationExactForLinear(self):
        grid, box = grid_bridge.SampleFieldsOnGrid(
            self.model_part, [("PRESSURE", "node_historical")], (4, 4, 4))
        rng = numpy.random.default_rng(0)
        points = rng.uniform(0.05, 0.95, size=(20, 3))
        values = grid_bridge.InterpolateGridAtPoints(grid, box, points)
        self.assertTrue(numpy.allclose(values[:, 0], _ExactLinear(points), atol=1e-12))

    def test_ScatterRoundTripExactForLinear(self):
        grid, box = grid_bridge.SampleFieldsOnGrid(
            self.model_part, [("PRESSURE", "node_historical")], (4, 4, 4))
        target = CreateStructuredTetModelPart(self.model, "Target", divisions=2,
                                              historical_variables=(Kratos.TEMPERATURE,))
        grid_bridge.ScatterGridToNodes(grid, box, target, [("TEMPERATURE", "node_historical")])
        for node in target.Nodes:
            expected = 1.0 + 2.0 * node.X + 3.0 * node.Y - node.Z
            self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.TEMPERATURE), expected, places=10)

    def test_MultiFieldChannelLayout(self):
        grid, _ = grid_bridge.SampleFieldsOnGrid(
            self.model_part, [("PRESSURE", "node_historical"), ("PRESSURE", "node_historical")], (3, 3, 3))
        self.assertEqual(grid.shape[0], 2)
        self.assertTrue(numpy.allclose(grid[0], grid[1]))

    def test_HexahedralFallbackIsExactOnUnsortedIds(self):
        # the per-point locator path (non-simplex geometries) - untested
        # before this - with node ids created out of order and a linear
        # field reproduced exactly. Note what this CANNOT discriminate:
        # Kratos keeps ModelPart.Nodes id-sorted whatever the creation
        # order, so the lookup's argsort is the identity on every real
        # part (a mutation dropping it survives here; TestRowsOfIds pins
        # it on an unsorted array).
        model = Kratos.Model()
        part = model.CreateModelPart("Hexes")
        part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        properties = part.CreateNewProperties(1)
        corners = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
                   (2, 0, 0), (2, 1, 0), (2, 0, 1), (2, 1, 1)]
        ids = [50, 10, 70, 20, 90, 30, 80, 40, 60, 100, 15, 25]  # deliberately unsorted
        for node_id, (x, y, z) in zip(ids, corners):
            node = part.CreateNewNode(node_id, float(x), float(y), float(z))
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0 + 2.0 * x - 3.0 * y + 0.5 * z)
        part.CreateNewElement("Element3D8N", 1, ids[:8], properties)
        part.CreateNewElement("Element3D8N", 2, [ids[1], ids[8], ids[9], ids[2], ids[5], ids[10], ids[11], ids[6]],
                              properties)
        grid, box = grid_bridge.SampleFieldsOnGrid(
            part, [("PRESSURE", "node_historical")], (5, 3, 3),
            bounding_box=(numpy.array([0.1, 0.1, 0.1]), numpy.array([1.9, 0.9, 0.9])))
        points = grid_bridge._GridPointCoordinates((5, 3, 3), box)
        expected = 1.0 + 2.0 * points[:, 0] - 3.0 * points[:, 1] + 0.5 * points[:, 2]
        numpy.testing.assert_allclose(grid.reshape(-1), expected, rtol=1e-10)

    def test_GaussPointLocationRejected(self):
        with self.assertRaisesRegex(ValueError, "nodal locations only"):
            grid_bridge.SampleFieldsOnGrid(
                self.model_part, [("PRESSURE", "element_gauss_point")], (3, 3, 3))

    def test_InvalidGridShapeRejected(self):
        with self.assertRaisesRegex(ValueError, "grid_shape"):
            grid_bridge.SampleFieldsOnGrid(
                self.model_part, [("PRESSURE", "node_historical")], (1, 4, 4))


try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import torch
    from physicsnemo.nn.functional import derivatives as _pn_derivatives  # noqa: F401
    have_grid_derivatives = True
except ImportError:
    have_grid_derivatives = False


@KratosUnittest.skipUnless(have_grid_derivatives, "Missing required python modules: torch, physicsnemo.")
class TestGridDerivatives(KratosUnittest.TestCase):
    """physicsnemo grid-derivative operators through ComputeGridDerivatives."""

    def _PeriodicField(self, n=64):
        axis = torch.arange(n, dtype=torch.float64) * (2.0 * torch.pi / n)
        x, y = torch.meshgrid(axis, axis, indexing="ij")
        return torch.sin(x) * torch.cos(y), x, y, float(axis[1])

    def test_UniformPeriodicAccuracy(self):
        field, x, y, spacing = self._PeriodicField()
        grid = field[None]  # (1, H, W)
        result, interior = grid_bridge.ComputeGridDerivatives(grid, Kratos.Parameters(
            '{"operator": "uniform", "spacing": [%.17g], "derivative_orders": [1]}' % spacing))
        self.assertEqual(tuple(result.shape), (1, 2, 64, 64))
        self.assertEqual(interior, (slice(None), slice(None)))
        self.assertLess(float((result[0, 0] - torch.cos(x) * torch.cos(y)).abs().max()), 2e-3)
        self.assertLess(float((result[0, 1] + torch.sin(x) * torch.sin(y)).abs().max()), 2e-3)
        # order-4 stencil: ~1e-5
        result4, _ = grid_bridge.ComputeGridDerivatives(grid, Kratos.Parameters(
            '{"operator": "uniform", "spacing": [%.17g], "order": 4}' % spacing))
        self.assertLess(float((result4[0, 0] - torch.cos(x) * torch.cos(y)).abs().max()), 1e-4)

    def test_SpectralPeriodicAccuracy(self):
        field, x, y, _ = self._PeriodicField()
        result, _ = grid_bridge.ComputeGridDerivatives(field[None], Kratos.Parameters(
            '{"operator": "spectral", "lengths": [%.17g]}' % (2.0 * torch.pi)))
        self.assertLess(float((result[0, 0] - torch.cos(x) * torch.cos(y)).abs().max()), 1e-4)

    def test_DerivativeChannelOrderingPin(self):
        field, x, y, spacing = self._PeriodicField()
        result, _ = grid_bridge.ComputeGridDerivatives(field[None], Kratos.Parameters(
            '{"operator": "spectral", "lengths": [%.17g],'
            ' "derivative_orders": [1, 2], "include_mixed": true}' % (2.0 * torch.pi)))
        # pinned order: dx, dy, dxx, dyy, dxy
        self.assertEqual(tuple(result.shape), (1, 5, 64, 64))
        sin, cos = torch.sin, torch.cos
        expected = [cos(x) * cos(y), -sin(x) * sin(y),
                    -sin(x) * cos(y), -sin(x) * cos(y), -cos(x) * sin(y)]
        for channel, reference in enumerate(expected):
            self.assertLess(float((result[0, channel] - reference).abs().max()), 1e-3,
                            msg=f"derivative channel {channel}")

    def test_TrimMakesNonPeriodicInteriorExact(self):
        n = 17
        axis = torch.linspace(0.0, 1.0, n, dtype=torch.float64)
        x, y = torch.meshgrid(axis, axis, indexing="ij")
        field = (x**2 + 3.0 * y)[None]  # non-periodic polynomial
        spacing = float(axis[1])
        result, interior = grid_bridge.ComputeGridDerivatives(field, Kratos.Parameters(
            '{"operator": "uniform", "spacing": [%.17g], "boundary": "trim"}' % spacing))
        self.assertEqual(tuple(result.shape), (1, 2, n - 2, n - 2))
        self.assertEqual(interior, (slice(1, -1), slice(1, -1)))
        self.assertLess(float((result[0, 0] - (2.0 * x)[interior]).abs().max()), 1e-10)
        self.assertLess(float((result[0, 1] - 3.0).abs().max()), 1e-10)
        # without trim the wrapped boundary layer is garbage
        untrimmed, _ = grid_bridge.ComputeGridDerivatives(field, Kratos.Parameters(
            '{"operator": "uniform", "spacing": [%.17g]}' % spacing))
        self.assertGreater(float((untrimmed[0, 0] - 2.0 * x).abs().max()), 1.0)

    def test_RectilinearNonUniformAxis(self):
        n = 33
        axis = torch.linspace(0.0, 1.0, n, dtype=torch.float64)**2
        field = (3.0 * axis)[None]  # f = 3x on a stretched axis
        settings = Kratos.Parameters('{"operator": "rectilinear", "boundary": "trim", "coordinates": []}')
        settings["coordinates"].Append(Kratos.Vector(axis.tolist()))
        result, interior = grid_bridge.ComputeGridDerivatives(field, settings)
        # upstream's nonuniform stencil accumulates in float32 internally
        self.assertLess(float((result[0, 0] - 3.0).abs().max()), 1e-5)

    def test_SpectralRefusesTrim(self):
        field = torch.zeros(1, 8, 8, dtype=torch.float64)
        with self.assertRaisesRegex(ValueError, "Spectral"):
            grid_bridge.ComputeGridDerivatives(field, Kratos.Parameters(
                '{"operator": "spectral", "boundary": "trim"}'))

    def test_Differentiable(self):
        field, _, _, spacing = self._PeriodicField(16)
        grid = field[None].clone().requires_grad_(True)
        result, _ = grid_bridge.ComputeGridDerivatives(grid, Kratos.Parameters(
            '{"operator": "uniform", "spacing": [%.17g]}' % spacing))
        result.square().sum().backward()
        self.assertIsNotNone(grid.grad)
        self.assertGreater(float(grid.grad.abs().sum()), 0.0)


class TestFlatBoundingBoxInterpolation(KratosUnittest.TestCase):
    """A planar Kratos case has every node at the same Z, so its bounding
    box is FLAT in that axis. Dividing by a zero extent gives NaN, and NaN
    cast to an index is INT_MIN - which used to crash the interpolation
    with an out-of-bounds index instead of returning the field."""

    def setUp(self):
        self.bounding_box = (numpy.array([0.0, 0.0, 0.0]), numpy.array([1.0, 1.0, 0.0]))
        # the two Z slices of a flat box are sampled at the same place, so
        # they carry identical values
        plane = numpy.array([[1.0, 2.0], [3.0, 4.0]])
        self.grid = numpy.stack([plane, plane], axis=-1)[None]  # (1, 2, 2, 2)

    def test_NumpyPathInterpolatesOnAFlatAxis(self):
        values = grid_bridge.InterpolateGridAtPoints(
            self.grid, self.bounding_box,
            numpy.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.5, 0.0, 0.0]]))
        numpy.testing.assert_allclose(values.ravel(), [1.0, 4.0, 2.0], atol=1e-12)

    def test_PointsOffTheFlatPlaneStillLand(self):
        """Nothing constrains a query to sit exactly on the plane; it must
        clamp onto it rather than produce NaN."""
        values = grid_bridge.InterpolateGridAtPoints(
            self.grid, self.bounding_box, numpy.array([[0.0, 0.0, 0.7]]))
        numpy.testing.assert_allclose(values.ravel(), [1.0], atol=1e-12)


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestTorchInterpolation(KratosUnittest.TestCase):
    def test_MatchesTheNumpyPathAndFlatAxesToo(self):
        rng = numpy.random.default_rng(0)
        grid = rng.standard_normal((2, 4, 5, 6))
        for bounding_box in ((numpy.array([0.0, 0.0, 0.0]), numpy.array([1.0, 2.0, 3.0])),
                             (numpy.array([0.0, 0.0, 0.0]), numpy.array([1.0, 2.0, 0.0]))):
            points = rng.uniform(0.0, 1.0, size=(20, 3)) * numpy.maximum(
                bounding_box[1] - bounding_box[0], 1e-9)
            reference = grid_bridge.InterpolateGridAtPoints(grid, bounding_box, points)
            values = grid_bridge.InterpolateGridAtPointsTorch(
                torch.tensor(grid), bounding_box, points)
            numpy.testing.assert_allclose(values.numpy(), reference, atol=1e-12)

    def test_GradientsReachTheGrid(self):
        grid = torch.arange(24.0, dtype=torch.float64).reshape(1, 2, 3, 4).requires_grad_(True)
        bounding_box = (numpy.array([0.0, 0.0, 0.0]), numpy.array([1.0, 1.0, 1.0]))
        points = numpy.array([[0.3, 0.4, 0.5], [0.9, 0.1, 0.2]])
        grid_bridge.InterpolateGridAtPointsTorch(grid, bounding_box, points).sum().backward()
        # trilinear weights sum to one per query point
        self.assertAlmostEqual(float(grid.grad.sum()), float(len(points)), places=10)


if __name__ == '__main__':
    KratosUnittest.main()
