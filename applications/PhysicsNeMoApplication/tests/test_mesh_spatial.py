"""Tests for signed distance fields on Kratos meshes and the grid
divergence/curl/laplacian operators (both physicsnemo >= 2.2)."""

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import spatial

from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    from physicsnemo.mesh.spatial import signed_distance_field  # noqa: F401
    have_sdf = True
except ImportError:
    have_sdf = False

try:
    import torch
    from physicsnemo.nn.functional import derivatives as _pn_derivatives
    have_vector_ops = hasattr(_pn_derivatives, "uniform_grid_divergence")
except ImportError:
    have_vector_ops = False

_MISSING = "Missing required python modules: torch, physicsnemo >= 2.2."


@KratosUnittest.skipUnless(have_sdf, _MISSING)
class TestSignedDistance(KratosUnittest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Sdf", 4, historical_variables=(Kratos.PRESSURE, Kratos.DISTANCE))

    def test_SphereDistancesMatchTheAnalyticValue(self):
        from physicsnemo.mesh.primitives.surfaces import sphere_icosahedral
        sphere = sphere_icosahedral.load(subdivisions=4)
        queries = numpy.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0],
                               [2.0, 0.0, 0.0], [0.0, 3.0, 0.0]])
        distance = spatial.ComputeSignedDistance(sphere, queries)

        # analytic: |x| - R with R = 1, negative inside
        expected = numpy.linalg.norm(queries, axis=1) - 1.0
        numpy.testing.assert_allclose(distance, expected, atol=2e-3)
        self.assertLess(distance[0], 0.0)   # centre is inside
        self.assertGreater(distance[2], 0.0)  # outside
        self.assertEqual(distance.dtype, numpy.float64)

    def test_VolumeMeshIsReducedToItsBoundary(self):
        mesh, _ = domain_mesh_builder.BuildMesh(self.model_part)
        # a tetrahedral mesh is codimension 0 and upstream rejects it outright
        self.assertEqual(int(mesh.cells.shape[1]) - 1, 3)
        surface = spatial.BoundarySurface(mesh)
        self.assertEqual(int(surface.cells.shape[1]) - 1, 2)
        # ... and passing a surface through is a no-op
        self.assertIs(spatial.BoundarySurface(surface), surface)

    def test_ExtractedBoundaryIsOrientedOutward(self):
        """The sign of the whole field depends on this: upstream's own
        boundary extraction winds triangles inconsistently (signed volume 0),
        which makes interior points read as outside."""
        mesh, _ = domain_mesh_builder.BuildMesh(self.model_part)
        surface = spatial.BoundarySurface(mesh)
        points = torch.as_tensor(surface.points).double()
        cells = torch.as_tensor(surface.cells).long()
        a, b, c = points[cells[:, 0]], points[cells[:, 1]], points[cells[:, 2]]
        signed_volume = float(
            torch.einsum("ij,ij->i", a, torch.cross(b, c, dim=1)).sum() / 6.0)
        self.assertAlmostEqual(signed_volume, 1.0, places=9)   # the unit cube

        unoriented = mesh.get_boundary_mesh()
        points = torch.as_tensor(unoriented.points).double()
        cells = torch.as_tensor(unoriented.cells).long()
        a, b, c = points[cells[:, 0]], points[cells[:, 1]], points[cells[:, 2]]
        self.assertLess(abs(float(
            torch.einsum("ij,ij->i", a, torch.cross(b, c, dim=1)).sum() / 6.0)), 1e-9)

    def test_NodalDistanceIsZeroOnTheBoundaryAndNegativeInside(self):
        distance = spatial.ComputeNodalSignedDistance(self.model_part)
        coordinates = numpy.array([[n.X, n.Y, n.Z] for n in self.model_part.Nodes])
        on_boundary = numpy.any(
            (coordinates <= 1e-12) | (coordinates >= 1.0 - 1e-12), axis=1)

        self.assertLess(numpy.abs(distance[on_boundary]).max(), 1e-6)
        self.assertTrue((distance[~on_boundary] < 0.0).all())
        # the deepest interior node of a unit cube is at distance 0.5
        self.assertAlmostEqual(float(distance.min()), -0.5, places=3)

    def test_WriteSignedDistanceFieldFeedsTheNormalGathers(self):
        written = spatial.WriteSignedDistanceField(self.model_part, Kratos.Parameters("""{
            "output_variable" : "DISTANCE",
            "output_location" : "node_non_historical"
        }"""))
        stored = numpy.array([n.GetValue(Kratos.DISTANCE) for n in self.model_part.Nodes])
        numpy.testing.assert_allclose(stored, written, atol=1e-12)

        # the point of writing a variable: every existing gather picks it up
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
        node_features, _, _, _ = graph_bridge.BuildGraph(
            self.model_part, (("DISTANCE", "node_non_historical"),))
        numpy.testing.assert_allclose(node_features[:, 0], written, atol=1e-12)

    def test_NarrowBandReturnsNaNOutsideIt(self):
        from physicsnemo.mesh.primitives.surfaces import sphere_icosahedral
        sphere = sphere_icosahedral.load(subdivisions=3)
        queries = numpy.array([[1.05, 0.0, 0.0], [8.0, 0.0, 0.0]])
        distance = spatial.ComputeSignedDistance(sphere, queries, max_dist=0.5)
        self.assertTrue(numpy.isfinite(distance[0]))
        self.assertTrue(numpy.isnan(distance[1]))

    def test_SdfSampledOnTheGrid(self):
        grid, bounding_box = spatial.SampleSignedDistanceOnGrid(self.model_part, (6, 6, 6))
        self.assertEqual(grid.shape, (1, 6, 6, 6))
        # corners of the lattice sit on the cube's surface, the centre is deepest
        self.assertAlmostEqual(float(grid[0, 0, 0, 0]), 0.0, places=5)
        self.assertLess(float(grid[0, 2, 2, 2]), -0.2)
        # it concatenates straight onto a sampled field grid
        fields, _ = grid_bridge.SampleFieldsOnGrid(
            self.model_part, [("PRESSURE", "node_historical")], (6, 6, 6),
            bounding_box=bounding_box)
        combined = numpy.concatenate([fields, grid], axis=0)
        self.assertEqual(combined.shape, (2, 6, 6, 6))


@KratosUnittest.skipUnless(have_vector_ops, _MISSING)
class TestGridVectorOperators(KratosUnittest.TestCase):

    def _PeriodicAxes(self, n=48):
        axis = torch.arange(n, dtype=torch.float64) * (2.0 * torch.pi / n)
        x, y = torch.meshgrid(axis, axis, indexing="ij")
        return x, y, float(axis[1])

    def test_DivergenceOfADivergenceFreeField(self):
        x, y, spacing = self._PeriodicAxes()
        field = torch.stack([torch.sin(x) * torch.cos(y),
                             -torch.cos(x) * torch.sin(y)], dim=0)
        result, interior = grid_bridge.ComputeGridVectorOperator(field, Kratos.Parameters(
            '{"operation": "divergence", "spacing": [%.17g]}' % spacing))
        self.assertEqual(tuple(result.shape), (48, 48))
        self.assertLess(float(result.abs().max()), 1e-3)

    def test_DivergenceOfALinearFieldIsExact(self):
        n = 24
        axis = torch.linspace(0.0, 1.0, n, dtype=torch.float64)
        x, y = torch.meshgrid(axis, axis, indexing="ij")
        field = torch.stack([2.0 * x, 3.0 * y], dim=0)   # divergence = 5
        result, interior = grid_bridge.ComputeGridVectorOperator(field, Kratos.Parameters(
            '{"operation": "divergence", "spacing": [%.17g], "boundary": "trim"}'
            % float(axis[1])))
        self.assertEqual(tuple(result.shape), (n - 2, n - 2))
        self.assertLess(float((result - 5.0).abs().max()), 1e-9)

    def test_CurlOfARigidRotationIsTwiceTheRate(self):
        n = 24
        axis = torch.linspace(0.0, 1.0, n, dtype=torch.float64)
        x, y = torch.meshgrid(axis, axis, indexing="ij")
        field = torch.stack([-y, x], dim=0)   # omega = 1 -> curl = 2
        result, _ = grid_bridge.ComputeGridVectorOperator(field, Kratos.Parameters(
            '{"operation": "curl", "spacing": [%.17g], "boundary": "trim"}' % float(axis[1])))
        self.assertEqual(tuple(result.shape), (n - 2, n - 2))   # 2D curl is scalar
        self.assertLess(float((result - 2.0).abs().max()), 1e-9)

    def test_Curl3DKeepsItsVectorShape(self):
        n = 12
        axis = torch.linspace(0.0, 1.0, n, dtype=torch.float64)
        x, y, z = torch.meshgrid(axis, axis, axis, indexing="ij")
        field = torch.stack([-y, x, torch.zeros_like(z)], dim=0)
        result, _ = grid_bridge.ComputeGridVectorOperator(field, Kratos.Parameters(
            '{"operation": "curl", "spacing": [%.17g], "boundary": "trim"}' % float(axis[1])))
        self.assertEqual(tuple(result.shape), (3, n - 2, n - 2, n - 2))
        self.assertLess(float((result[2] - 2.0).abs().max()), 1e-9)

    def test_LaplacianOfAHarmonicFieldVanishes(self):
        n = 24
        axis = torch.linspace(0.0, 1.0, n, dtype=torch.float64)
        x, y = torch.meshgrid(axis, axis, indexing="ij")
        harmonic = x * x - y * y    # laplacian = 0
        result, _ = grid_bridge.ComputeGridVectorOperator(harmonic, Kratos.Parameters(
            '{"operation": "laplacian", "spacing": [%.17g], "boundary": "trim"}'
            % float(axis[1])))
        self.assertEqual(tuple(result.shape), (n - 2, n - 2))
        self.assertLess(float(result.abs().max()), 1e-8)

    def test_LaplacianAppliesPerChannel(self):
        n = 20
        axis = torch.linspace(0.0, 1.0, n, dtype=torch.float64)
        x, y = torch.meshgrid(axis, axis, indexing="ij")
        stacked = torch.stack([x * x, y * y], dim=0)   # laplacians 2 and 2
        result, _ = grid_bridge.ComputeGridVectorOperator(stacked, Kratos.Parameters(
            '{"operation": "laplacian", "spacing": [%.17g], "boundary": "trim",'
            ' "has_channel_axis": true}' % float(axis[1])))
        self.assertEqual(tuple(result.shape), (2, n - 2, n - 2))
        self.assertLess(float((result - 2.0).abs().max()), 1e-8)

    def test_Validation(self):
        field = torch.zeros(2, 8, 8, dtype=torch.float64)
        with self.assertRaisesRegex(ValueError, "Unknown operation"):
            grid_bridge.ComputeGridVectorOperator(
                field, Kratos.Parameters('{"operation": "gradient"}'))
        with self.assertRaisesRegex(ValueError, "spectral"):
            grid_bridge.ComputeGridVectorOperator(
                field, Kratos.Parameters('{"operator": "spectral"}'))
        # a channel count that does not match the spatial rank is the classic
        # mistake, since gradients take a scalar field but these take a vector
        mismatched = torch.zeros(3, 8, 8, dtype=torch.float64)
        with self.assertRaisesRegex(ValueError, "channel count"):
            grid_bridge.ComputeGridVectorOperator(
                mismatched, Kratos.Parameters('{"operation": "divergence"}'))


if __name__ == '__main__':
    KratosUnittest.main()
