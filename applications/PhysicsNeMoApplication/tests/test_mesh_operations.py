"""Tests for the physicsnemo.mesh operations the bridge now exposes:
curvature, repair, subdivision, smoothing, extrusion, moments, the two
mesh-level deformers, and torch-native sampling over a BVH."""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import physicsnemo.mesh.curvature  # noqa: F401
    import physicsnemo.mesh.sampling  # noqa: F401
    have_mesh_operations = True
except ImportError:
    have_mesh_operations = False


def _UnitSphere(levels=3):
    """A subdivided octahedron projected onto the unit sphere: a surface
    whose curvature and area are known in closed form."""
    from physicsnemo.mesh import Mesh
    from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import operations

    points = torch.tensor([[1.0, 0, 0], [-1.0, 0, 0], [0, 1.0, 0],
                           [0, -1.0, 0], [0, 0, 1.0], [0, 0, -1.0]])
    cells = torch.tensor([[0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4],
                          [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5]])
    mesh = Mesh(points=points, cells=cells)
    for _ in range(levels):
        mesh = operations.SubdivideMesh(mesh, "loop")
        mesh = Mesh(points=mesh.points / mesh.points.norm(dim=1, keepdim=True),
                    cells=mesh.cells)
    return mesh


@KratosUnittest.skipUnless(have_torch and have_mesh_operations,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestMeshOperations(KratosUnittest.TestCase):
    def test_CurvatureOfTheUnitSphere(self):
        """Both curvatures of a unit sphere are 1, which is what makes this
        the test that the bridge hands the surface over correctly."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import operations

        sphere = _UnitSphere()
        mean = operations.ComputeCurvature(sphere, "mean")
        gaussian = operations.ComputeCurvature(sphere, "gaussian")
        self.assertAlmostEqual(float(mean.median()), 1.0, places=2)
        self.assertAlmostEqual(float(gaussian.median()), 1.0, places=1)

        with self.assertRaisesRegex(ValueError, "curvature"):
            operations.ComputeCurvature(sphere, "principal")

    def test_SubdivisionQuadruplesTheCells(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import operations

        sphere = _UnitSphere(levels=1)
        for scheme in ("loop", "butterfly", "linear"):
            refined = operations.SubdivideMesh(sphere, scheme)
            self.assertEqual(refined.cells.shape[0], 4 * sphere.cells.shape[0])
        twice = operations.SubdivideMesh(sphere, "linear", levels=2)
        self.assertEqual(twice.cells.shape[0], 16 * sphere.cells.shape[0])

        with self.assertRaisesRegex(ValueError, "subdivision scheme"):
            operations.SubdivideMesh(sphere, "catmull_clark")
        with self.assertRaisesRegex(ValueError, "levels"):
            operations.SubdivideMesh(sphere, "loop", levels=0)

    def test_LinearSubdivisionKeepsThePointsPutAndLoopDoesNot(self):
        """The distinction that matters when refining a mesh you intend to
        keep: a smoothing subdivision MOVES the surface."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import operations

        sphere = _UnitSphere(levels=2)
        linear = operations.SubdivideMesh(sphere, "linear")
        loop = operations.SubdivideMesh(sphere, "loop")
        original = sphere.points.shape[0]
        torch.testing.assert_close(linear.points[:original], sphere.points)
        self.assertGreater(
            float((loop.points[:original] - sphere.points).abs().max()), 1e-6)

    def test_IntegrateMomentOfOnesIsTheArea(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import operations

        sphere = _UnitSphere()
        ones = torch.ones(sphere.cells.shape[0])
        area = float(operations.IntegrateMoment(sphere, ones, ones))
        self.assertAlmostEqual(area, 4.0 * numpy.pi, delta=0.2)  # a tessellation undershoots
        self.assertLess(area, 4.0 * numpy.pi)

    def test_SmoothingRelaxesANoisySphere(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import operations
        from physicsnemo.mesh import Mesh

        sphere = _UnitSphere()
        torch.manual_seed(0)
        noisy = Mesh(points=sphere.points + 0.05 * torch.randn_like(sphere.points),
                     cells=sphere.cells)
        smoothed = operations.SmoothMesh(noisy, n_iter=20, relaxation_factor=0.2)
        spread = lambda mesh: float(mesh.points.norm(dim=1).std())
        self.assertLess(spread(smoothed), spread(noisy))

    def test_RepairReportsWhatItChanged(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import operations
        from physicsnemo.mesh import Mesh

        sphere = _UnitSphere(levels=1)
        # a duplicated point and a degenerate cell referencing it
        points = torch.cat([sphere.points, sphere.points[:1]], dim=0)
        duplicate = points.shape[0] - 1
        cells = torch.cat(
            [sphere.cells, torch.tensor([[0, 0, duplicate]])], dim=0)
        broken = Mesh(points=points, cells=cells)

        repaired, report = operations.RepairMesh(broken)
        self.assertIn("degenerates", report)
        self.assertLess(repaired.cells.shape[0], broken.cells.shape[0])


@KratosUnittest.skipUnless(have_torch and have_mesh_operations,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestCurvatureOnAModelPart(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=3,
            historical_variables=(Kratos.PRESSURE,))

    def test_CurvatureLandsOnTheNodes(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import operations

        operations.WriteCurvatureFields(self.model_part, Kratos.Parameters("""{
            "mean_curvature_variable"     : "NODAL_PAUX",
            "gaussian_curvature_variable" : "NODAL_AREA"
        }"""))
        mean = numpy.array([node.GetValue(Kratos.NODAL_PAUX)
                            for node in self.model_part.Nodes])
        gaussian = numpy.array([node.GetValue(Kratos.NODAL_AREA)
                                for node in self.model_part.Nodes])
        self.assertEqual(len(mean), self.model_part.NumberOfNodes())
        # upstream returns NaN at a node on no boundary facet - every
        # interior node of a volume mesh - because its vertex area is zero
        self.assertTrue(numpy.isfinite(mean).all())
        self.assertTrue(numpy.isfinite(gaussian).all())
        interior = numpy.array([
            1e-12 < node.X < 1.0 - 1e-12 and 1e-12 < node.Y < 1.0 - 1e-12
            and 1e-12 < node.Z < 1.0 - 1e-12 for node in self.model_part.Nodes])
        self.assertGreater(interior.sum(), 0)
        numpy.testing.assert_array_equal(mean[interior], 0.0)

        # a cube is flat on its faces and sharp at its corners: the Gaussian
        # curvature concentrates on the eight corner nodes
        corners = numpy.array([
            abs(node.X - round(node.X)) < 1e-12 and abs(node.Y - round(node.Y)) < 1e-12
            and abs(node.Z - round(node.Z)) < 1e-12
            and node.X in (0.0, 1.0) and node.Y in (0.0, 1.0) and node.Z in (0.0, 1.0)
            for node in self.model_part.Nodes])
        self.assertEqual(corners.sum(), 8)
        self.assertGreater(numpy.abs(gaussian[corners]).mean(),
                           numpy.abs(gaussian[~corners]).mean())

    def test_WritingNothingRaises(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import operations

        with self.assertRaisesRegex(ValueError, "nothing to write"):
            operations.WriteCurvatureFields(self.model_part, Kratos.Parameters("{}"))


@KratosUnittest.skipUnless(have_torch and have_mesh_operations,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestExtrusion(KratosUnittest.TestCase):
    """A 2-D Kratos case swept into a 3-D one, as real entities."""

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Planar")
        self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
        properties = self.model_part.CreateNewProperties(1)
        for node_id, (x, y) in enumerate(
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], start=1):
            self.model_part.CreateNewNode(node_id, x, y, 0.0)
        self.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)
        self.model_part.CreateNewElement("Element2D3N", 2, [1, 3, 4], properties)

    def test_ASweptSquareIsAVolume(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import operations

        swept = operations.ExtrudeModelPart(
            self.model, self.model_part, [0.0, 0.0, 2.0], "Swept")
        self.assertEqual(swept.NumberOfNodes(), 8)     # 4 points, two layers
        # upstream decomposes the sweep into SIMPLICES, not prisms: each
        # swept triangle arrives as three tetrahedra
        self.assertEqual(swept.NumberOfElements(), 6)
        z_values = sorted({round(node.Z, 10) for node in swept.Nodes})
        self.assertEqual(z_values, [0.0, 2.0])

    def test_AnInPlaneVectorIsRefused(self):
        """Sweeping along the source's own plane sweeps nothing, and
        upstream returns the flat result rather than complaining - so the
        measure is what says so."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import operations

        with self.assertRaisesRegex(ValueError, "plane"):
            operations.ExtrudeModelPart(
                self.model, self.model_part, [1.0, 0.0, 0.0], "Flat")

    def test_TheSweptVolumeIsAreaTimesLength(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
            domain_mesh_builder, operations)

        swept = operations.ExtrudeModelPart(
            self.model, self.model_part, [0.0, 0.0, 2.0], "Measured")
        mesh, _ = domain_mesh_builder.BuildMesh(swept)
        ones = torch.ones(mesh.cells.shape[0], dtype=mesh.points.dtype)
        self.assertAlmostEqual(
            float(operations.IntegrateMoment(mesh, ones, ones)), 2.0, places=8)


@KratosUnittest.skipUnless(have_torch and have_mesh_operations,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestMeshDeformers(KratosUnittest.TestCase):
    """shrinkwrap and sobolev_deform, which take a MESH rather than points."""

    def test_ShrinkwrapProjectsOntoTheTarget(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import deformation
        from physicsnemo.mesh import Mesh

        target = _UnitSphere(levels=3)
        source = _UnitSphere(levels=2)
        scattered = Mesh(points=source.points * 1.4, cells=source.cells)
        self.assertGreater(float(scattered.points.norm(dim=1).mean()), 1.3)

        wrapped = deformation.DeformMesh(scattered, "shrinkwrap", target=target)
        radii = wrapped.points.norm(dim=1)
        self.assertAlmostEqual(float(radii.mean()), 1.0, places=2)
        self.assertLess(float(radii.std()), 0.02)

    def test_SobolevSpreadsADisplacementAndPinsFixedPoints(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import deformation

        sphere = _UnitSphere(levels=3)
        # push a single vertex outward and let the filter spread it
        displacement = torch.zeros_like(sphere.points)
        displacement[0, 0] = 1.0
        deformed = deformation.DeformMesh(
            sphere, "sobolev", displacement=displacement, length_scale=0.5)
        moved = (deformed.points - sphere.points).norm(dim=1)
        self.assertGreater(float(moved[0]), 0.0)
        # the neighbourhood moved too - that is the whole point of the filter
        self.assertGreater(int((moved > 1e-6).sum()), 1)

    def test_MethodValidation(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import deformation

        sphere = _UnitSphere(levels=1)
        with self.assertRaisesRegex(ValueError, "mesh deformation"):
            deformation.DeformMesh(sphere, "ffd")
        with self.assertRaisesRegex(ValueError, "target"):
            deformation.DeformMesh(sphere, "shrinkwrap")
        with self.assertRaisesRegex(ValueError, "length_scale"):
            deformation.DeformMesh(
                sphere, "sobolev", displacement=torch.zeros_like(sphere.points))


@KratosUnittest.skipUnless(have_torch and have_mesh_operations,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestBvhSampling(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=4,
            historical_variables=(Kratos.PRESSURE,))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 2.0 * node.X + 3.0 * node.Y)

    def test_SamplingALinearFieldIsExact(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import sampling

        points = numpy.array([[0.25, 0.25, 0.5], [0.5, 0.75, 0.25], [0.9, 0.1, 0.9]])
        values, mesh = sampling.SampleModelPartAtPoints(
            self.model_part, [("PRESSURE", "node_historical")], points)
        expected = 2.0 * points[:, 0] + 3.0 * points[:, 1]
        numpy.testing.assert_allclose(
            values.detach().numpy().ravel(), expected, atol=1e-6)

        # the BVH is reusable across queries, which is the point of building it
        bvh = sampling.BuildBvh(mesh)
        again = sampling.SampleMeshAtPoints(mesh, points, bvh=bvh)
        numpy.testing.assert_allclose(
            numpy.asarray(again["PRESSURE"].detach().numpy()).ravel(), expected, atol=1e-6)

    def test_TheTwoGridBackendsAgree(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge

        specs = [("PRESSURE", "node_historical")]
        locator, box = grid_bridge.SampleFieldsOnGrid(
            self.model_part, specs, (6, 6, 6), backend="kratos")
        bvh, _ = grid_bridge.SampleFieldsOnGrid(
            self.model_part, specs, (6, 6, 6), bounding_box=box, backend="physicsnemo")
        self.assertEqual(locator.shape, bvh.shape)
        numpy.testing.assert_allclose(bvh, locator, atol=1e-6)

    def test_AnUnknownBackendRaises(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge

        with self.assertRaisesRegex(ValueError, "sampling backend"):
            grid_bridge.SampleFieldsOnGrid(
                self.model_part, [("PRESSURE", "node_historical")], (4, 4, 4),
                backend="octree")

    def test_SamplingCarriesGradientsToTheMeshData(self):
        """What the locator path cannot do at all."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
            domain_mesh_builder, sampling)

        mesh, _ = domain_mesh_builder.BuildMesh(self.model_part)
        mesh.point_data["u"] = torch.zeros(
            mesh.points.shape[0], dtype=torch.float64, requires_grad=True)
        sampled = sampling.SampleMeshAtPoints(
            mesh, numpy.array([[0.3, 0.4, 0.5]]), data_source="points")
        sampled["u"].sum().backward()
        self.assertIsNotNone(mesh.point_data["u"].grad)
        self.assertGreater(float(mesh.point_data["u"].grad.abs().sum()), 0.0)


if __name__ == '__main__':
    KratosUnittest.main()
