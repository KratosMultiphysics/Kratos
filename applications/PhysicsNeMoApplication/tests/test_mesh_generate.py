"""Tests for mesh generation and repair from implicit geometry, and for
materializing a generated mesh as real Kratos entities."""

import math

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import generate
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import spatial

from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    import physicsnemo.mesh
    from physicsnemo.mesh.generate import mesh_implicit_domain  # noqa: F401
    have_generate = True
except ImportError:
    have_generate = False

if have_generate:
    try:
        import tetgen  # noqa: F401
        have_tetgen = True
    except ImportError:
        have_tetgen = False
else:
    have_tetgen = False

have_meshing = kratos_utils.CheckIfApplicationsAvailable("MeshingApplication")
have_mmg = False
if have_meshing:
    try:
        import KratosMultiphysics.MeshingApplication as MeshingApplication
        have_mmg = hasattr(MeshingApplication, "MmgProcess3D")
    except ImportError:
        # Present in the install tree but not loadable (partial or
        # mismatched build); degrade to a skip instead of killing the module.
        have_meshing = False
have_mapping = kratos_utils.CheckIfApplicationsAvailable("MappingApplication")

_MISSING = "Missing required python modules: torch, physicsnemo >= 2.2."


def _SignedVolume(mesh):
    """Closed-surface signed volume; its sign reveals the triangle winding."""
    points = torch.as_tensor(mesh.points).double()
    cells = torch.as_tensor(mesh.cells).long()
    a, b, c = points[cells[:, 0]], points[cells[:, 1]], points[cells[:, 2]]
    return float(torch.einsum("ij,ij->i", a, torch.cross(b, c, dim=1)).sum() / 6.0)


def _TetVolume(mesh):
    points = torch.as_tensor(mesh.points).double()
    cells = torch.as_tensor(mesh.cells).long()
    a, b, c, d = (points[cells[:, i]] for i in range(4))
    return float((torch.einsum("ij,ij->i", b - a, torch.cross(c - a, d - a, dim=1))
                  / 6.0).abs().sum())


@KratosUnittest.skipUnless(have_generate, _MISSING)
class TestImplicitDomainGeneration(KratosUnittest.TestCase):

    def setUp(self):
        self.sphere = generate.SdfPrimitives()["sphere"]([0.0, 0.0, 0.0], 0.7)
        self.bounds = ([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])

    def test_GeneratedVolumeConvergesToTheAnalyticOne(self):
        exact = 4.0 / 3.0 * math.pi * 0.7 ** 3
        errors = []
        for h in (0.25, 0.15):
            mesh = generate.GenerateImplicitDomain(self.sphere, self.bounds, h)
            self.assertEqual(int(mesh.cells.shape[1]), 4)   # tets
            errors.append(abs(_TetVolume(mesh) - exact) / exact)
        self.assertLess(errors[1], errors[0])
        self.assertLess(errors[1], 0.05)

    def test_DiagnosticsReportAValidMesh(self):
        mesh, diagnostics = generate.GenerateImplicitDomain(
            self.sphere, self.bounds, 0.2, Kratos.Parameters('{"full_output": true}'))
        # q_min is NOT monotone in h, so the meaningful assertions are these
        self.assertTrue(diagnostics["all_volumes_positive"])
        self.assertTrue(diagnostics["boundary_closed_manifold"])
        self.assertGreater(diagnostics["q_median"], 0.5)
        self.assertLess(diagnostics["coverage_gap_h"], 1.5)

    def test_WorksInsideNoGrad(self):
        """Regression: upstream needs grad mode even though nothing trains.

        Called plainly inside torch.no_grad() the generator's boundary
        projection degenerates - it trips the coverage guard, or silently
        returns a worse mesh when the guard is disabled. Any deployment
        process may have wrapped the world in no_grad, so the wrapper must
        re-enable it.
        """
        with torch.no_grad():
            mesh, diagnostics = generate.GenerateImplicitDomain(
                self.sphere, self.bounds, 0.2, Kratos.Parameters('{"full_output": true}'))
        self.assertTrue(diagnostics["all_volumes_positive"])
        self.assertGreater(diagnostics["q_median"], 0.5)

        # ... and the bare upstream call in the same context does fail
        from physicsnemo.mesh.generate import mesh_implicit_domain
        with torch.no_grad():
            with self.assertRaises(ValueError):
                mesh_implicit_domain(self.sphere, self.bounds, 0.2, device="cpu")

    def test_DeterministicOnCpu(self):
        first = generate.GenerateImplicitDomain(self.sphere, self.bounds, 0.25)
        second = generate.GenerateImplicitDomain(self.sphere, self.bounds, 0.25)
        self.assertTrue(torch.equal(torch.as_tensor(first.points),
                                    torch.as_tensor(second.points)))

    def test_Validation(self):
        with self.assertRaisesRegex(ValueError, "reconnect"):
            generate.GenerateImplicitDomain(self.sphere, self.bounds, 0.25,
                                            Kratos.Parameters('{"reconnect": "bogus"}'))
        with self.assertRaisesRegex(ValueError, "high > low"):
            generate.GenerateImplicitDomain(
                self.sphere, ([1.0, 1.0, 1.0], [-1.0, -1.0, -1.0]), 0.25)


@KratosUnittest.skipUnless(have_generate, _MISSING)
class TestRefitToImplicit(KratosUnittest.TestCase):
    """The differentiable counterpart of GenerateImplicitDomain."""

    def setUp(self):
        self.primitives = generate.SdfPrimitives()
        self.bounds = ([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])

    def test_RefitMovesTheBoundaryOntoTheNewZeroSet(self):
        # Generate for one radius, then refit onto a slightly larger sphere:
        # the boundary must move toward the new zero set without the cell
        # count changing (topology is held fixed).
        mesh = generate.GenerateImplicitDomain(
            self.primitives["sphere"]([0.0, 0.0, 0.0], 0.6), self.bounds, 0.2)
        cells_before = int(mesh.cells.shape[0])
        volume_before = _TetVolume(mesh)

        target = self.primitives["sphere"]([0.0, 0.0, 0.0], 0.7)
        refitted = generate.RefitToImplicit(mesh, target, iters=3, bounds=self.bounds)

        self.assertEqual(int(refitted.cells.shape[0]), cells_before)
        self.assertGreater(_TetVolume(refitted), volume_before)

        exact = 4.0 / 3.0 * math.pi * 0.7 ** 3
        before = abs(volume_before - exact) / exact
        after = abs(_TetVolume(refitted) - exact) / exact
        self.assertLess(after, before)

    def test_RefitToItsOwnLevelSetBarelyMoves(self):
        # The near-identity case: refitting onto the SDF the mesh was
        # generated from should not drag it anywhere.
        sphere = self.primitives["sphere"]([0.0, 0.0, 0.0], 0.6)
        mesh = generate.GenerateImplicitDomain(sphere, self.bounds, 0.2)
        volume = _TetVolume(mesh)
        refitted = generate.RefitToImplicit(mesh, sphere, iters=2, bounds=self.bounds)
        self.assertAlmostEqual(_TetVolume(refitted), volume, delta=0.05 * volume)


@KratosUnittest.skipUnless(have_generate, _MISSING)
class TestSurfaceFromLevelSet(KratosUnittest.TestCase):

    def test_MarchingCubesIsOutwardWoundUnlikeBoundaryExtraction(self):
        n = 48
        axis = torch.linspace(-1.0, 1.0, n, dtype=torch.float64)
        x, y, z = torch.meshgrid(axis, axis, axis, indexing="ij")
        field = (x * x + y * y + z * z).sqrt() - 0.5

        surface = generate.SurfaceFromLevelSet(
            field, bounding_box=([-1.0] * 3, [1.0] * 3))
        exact = 4.0 / 3.0 * math.pi * 0.5 ** 3
        # a consistent OUTWARD winding gives +volume; an inconsistent one cancels
        self.assertAlmostEqual(_SignedVolume(surface), exact, delta=0.02 * exact)

        # the contrast that motivates using marching cubes for surfaces
        volume = generate.GenerateImplicitDomain(
            generate.SdfPrimitives()["sphere"]([0.0, 0.0, 0.0], 0.5),
            ([-1.0] * 3, [1.0] * 3), 0.2)
        self.assertAlmostEqual(_SignedVolume(volume.get_boundary_mesh()), 0.0, places=9)

    def test_AcceptsTheGridBridgeSdfLayout(self):
        model = Kratos.Model()
        model_part = CreateStructuredTetModelPart(
            model, "Sdf", 4, historical_variables=(Kratos.PRESSURE,))
        grid, bounding_box = spatial.SampleSignedDistanceOnGrid(model_part, (16, 16, 16))
        self.assertEqual(grid.shape[0], 1)     # (1, *grid) as shipped
        surface = generate.SurfaceFromLevelSet(grid, bounding_box=bounding_box)
        self.assertEqual(int(surface.cells.shape[1]), 3)
        self.assertGreater(int(surface.cells.shape[0]), 0)

    def test_Validation(self):
        with self.assertRaisesRegex(ValueError, "3D scalar field"):
            generate.SurfaceFromLevelSet(torch.zeros(4, 4, dtype=torch.float64))


@KratosUnittest.skipUnless(have_generate, _MISSING)
class TestFillBoundaryLoop(KratosUnittest.TestCase):

    def _CircleLoop(self, segments=32):
        angle = torch.linspace(0.0, 2.0 * math.pi, segments + 1)[:-1].double()
        points = torch.stack([angle.cos(), angle.sin()], dim=1)
        edges = torch.stack([torch.arange(segments),
                             (torch.arange(segments) + 1) % segments], dim=1).long()
        return physicsnemo.mesh.Mesh(points=points, cells=edges)

    def test_MinimumAngleGuaranteeIsMet(self):
        loop = self._CircleLoop()
        filled = generate.FillBoundaryLoop(loop)
        self.assertEqual(int(filled.cells.shape[1]), 3)

        points = torch.as_tensor(filled.points).double()
        cells = torch.as_tensor(filled.cells).long()
        smallest = math.inf
        for triangle in cells:
            a, b, c = (points[int(i)] for i in triangle)
            for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
                first, second = v - u, w - u
                cosine = float(first.dot(second)
                               / (first.norm() * second.norm()))
                smallest = min(smallest, math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
        self.assertGreaterEqual(smallest, 29.9)

        # the input loop survives as the leading rows
        self.assertTrue(torch.equal(points[:len(loop.points)],
                                    torch.as_tensor(loop.points).double()))

    def test_MaxCellSizeRefines(self):
        loop = self._CircleLoop()
        coarse = generate.FillBoundaryLoop(loop)
        fine = generate.FillBoundaryLoop(
            loop, Kratos.Parameters('{"max_cell_size": 0.005}'))
        self.assertGreater(int(fine.cells.shape[0]), int(coarse.cells.shape[0]))

    def test_ThreeDimensionalSurfaceRedirectedToTheTetrahedralFill(self):
        n = 24
        axis = torch.linspace(-1.0, 1.0, n, dtype=torch.float64)
        x, y, z = torch.meshgrid(axis, axis, axis, indexing="ij")
        surface = generate.SurfaceFromLevelSet((x * x + y * y + z * z).sqrt() - 0.5)
        with self.assertRaisesRegex(ValueError, "FillSurfaceWithTetrahedra"):
            generate.FillBoundaryLoop(surface)


@KratosUnittest.skipUnless(have_generate, _MISSING)
class TestPopulateModelPart(KratosUnittest.TestCase):

    def test_RoundTripsAKratosMesh(self):
        model = Kratos.Model()
        source = CreateStructuredTetModelPart(
            model, "Source", 3, historical_variables=(Kratos.PRESSURE,))
        mesh, _ = domain_mesh_builder.BuildMesh(source)

        rebuilt = generate.PopulateModelPartFromMesh(model, "Rebuilt", mesh)
        self.assertEqual(rebuilt.NumberOfNodes(), int(mesh.points.shape[0]))
        self.assertEqual(rebuilt.NumberOfElements(), int(mesh.cells.shape[0]))
        self.assertEqual(rebuilt.ProcessInfo[Kratos.DOMAIN_SIZE], 3)

        # geometry survives, and the 0-based -> 1-based shift is right
        rebuilt_mesh, _ = domain_mesh_builder.BuildMesh(rebuilt)
        numpy.testing.assert_allclose(
            numpy.asarray(rebuilt_mesh.points), numpy.asarray(mesh.points), atol=1e-12)
        numpy.testing.assert_array_equal(
            numpy.asarray(rebuilt_mesh.cells), numpy.asarray(mesh.cells))
        self.assertEqual(min(node.Id for node in rebuilt.Nodes), 1)

    def test_HistoricalVariablesAndBuffer(self):
        model = Kratos.Model()
        source = CreateStructuredTetModelPart(
            model, "Src2", 2, historical_variables=(Kratos.PRESSURE,))
        mesh, _ = domain_mesh_builder.BuildMesh(source)
        part = generate.PopulateModelPartFromMesh(model, "WithVars", mesh, Kratos.Parameters("""{
            "historical_variables" : ["TEMPERATURE"],
            "buffer_size"          : 2
        }"""))
        self.assertTrue(part.HasNodalSolutionStepVariable(Kratos.TEMPERATURE))
        self.assertEqual(part.GetBufferSize(), 2)
        for node in part.Nodes:
            node.SetSolutionStepValue(Kratos.TEMPERATURE, 1.5)
        self.assertAlmostEqual(
            next(iter(part.Nodes)).GetSolutionStepValue(Kratos.TEMPERATURE), 1.5)

    def test_Validation(self):
        model = Kratos.Model()
        source = CreateStructuredTetModelPart(
            model, "Src3", 2, historical_variables=(Kratos.PRESSURE,))
        mesh, _ = domain_mesh_builder.BuildMesh(source)
        generate.PopulateModelPartFromMesh(model, "Taken", mesh)
        with self.assertRaisesRegex(ValueError, "already exists"):
            generate.PopulateModelPartFromMesh(model, "Taken", mesh)
        with self.assertRaisesRegex(ValueError, "entity_type"):
            generate.PopulateModelPartFromMesh(
                model, "Bad", mesh, Kratos.Parameters('{"entity_type": "Nodes"}'))

    def test_GenerateModelPartEndToEnd(self):
        model = Kratos.Model()
        sphere = generate.SdfPrimitives()["sphere"]([0.0, 0.0, 0.0], 0.7)
        part = generate.GenerateModelPart(
            model, "Generated", sphere, ([-1.0] * 3, [1.0] * 3), 0.25)
        self.assertGreater(part.NumberOfNodes(), 0)
        self.assertGreater(part.NumberOfElements(), 0)
        self.assertEqual(part.ProcessInfo[Kratos.DOMAIN_SIZE], 3)
        radius = numpy.linalg.norm(
            numpy.array([[n.X, n.Y, n.Z] for n in part.Nodes]), axis=1)
        self.assertLess(radius.max(), 0.75)   # everything inside the sphere


@KratosUnittest.skipUnless(have_generate and have_mapping,
                           "Missing physicsnemo >= 2.2 or MappingApplication.")
class TestFieldTransferToGeneratedMesh(KratosUnittest.TestCase):

    def test_LinearFieldMapsOntoAGeneratedMesh(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import mapping_bridge
        model = Kratos.Model()
        source = CreateStructuredTetModelPart(
            model, "Old", 4, historical_variables=(Kratos.TEMPERATURE,))
        for node in source.Nodes:
            node.SetSolutionStepValue(Kratos.TEMPERATURE, 2.0 * node.X + 3.0 * node.Y)

        # a generated part covering the interior of the source cube
        # h must resolve the feature: 0.2 against a radius-0.3 sphere is only
        # ~3 cells across the diameter and trips the coverage guard
        sphere = generate.SdfPrimitives()["sphere"]([0.5, 0.5, 0.5], 0.3)
        destination = generate.GenerateModelPart(
            model, "New", sphere, ([0.1] * 3, [0.9] * 3), 0.1,
            populate_settings=Kratos.Parameters(
                '{"historical_variables": ["TEMPERATURE"], "buffer_size": 1}'))

        bridge = mapping_bridge.MappingBridge(
            source, destination, Kratos.Parameters('{"mapper_type": "nearest_element"}'))
        bridge.MapFields([("TEMPERATURE", "TEMPERATURE")])

        for node in destination.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE),
                2.0 * node.X + 3.0 * node.Y, places=8)


@KratosUnittest.skipUnless(have_generate and have_mmg,
                           "Missing physicsnemo >= 2.2 or MeshingApplication with MMG.")
class TestGeneratedMeshFeedsMmg(KratosUnittest.TestCase):

    def test_GeneratedMeshCanBeRemeshed(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import adaptive_remeshing
        model = Kratos.Model()
        sphere = generate.SdfPrimitives()["sphere"]([0.0, 0.0, 0.0], 0.7)
        part = generate.GenerateModelPart(
            model, "ForMmg", sphere, ([-1.0] * 3, [1.0] * 3), 0.25,
            populate_settings=Kratos.Parameters(
                '{"historical_variables": ["TEMPERATURE"], "buffer_size": 1}'))
        nodes_before = part.NumberOfNodes()

        sizes = adaptive_remeshing.ComputeTargetSizeField(
            part, numpy.full(nodes_before, 1e-3), Kratos.Parameters("""{
                "target_error": 1e-3, "exponent": 0.5,
                "minimal_size": 0.05, "maximal_size": 0.4
            }"""))
        adaptive_remeshing.RunMmgAdaptation(part, sizes)
        self.assertTrue(part.Is(Kratos.MODIFIED))
        self.assertGreater(part.NumberOfNodes(), 0)



def _CubeSurface():
    """Unit cube as 12 outward-wound triangles."""
    points = numpy.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                          [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
    triangles = numpy.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
                             [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
                             [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]])
    return points, triangles


def _LPrismSurface():
    """Non-convex L-shaped prism: base area 3, height 1, so volume 3."""
    base = numpy.array([[0, 0], [2, 0], [2, 1], [1, 1], [1, 2], [0, 2]], dtype=float)
    points = numpy.ascontiguousarray(numpy.vstack(
        [numpy.c_[base, numpy.zeros(6)], numpy.c_[base, numpy.ones(6)]]))
    cap = [(0, 1, 2), (0, 2, 3), (0, 3, 5), (3, 4, 5)]        # valid CCW triangulation
    triangles = [(c[0], c[2], c[1]) for c in cap]              # bottom, outward -z
    triangles += [(c[0] + 6, c[1] + 6, c[2] + 6) for c in cap]  # top, outward +z
    for i in range(6):
        j = (i + 1) % 6
        triangles += [(i, j, j + 6), (i, j + 6, i + 6)]
    return points, numpy.array(triangles)


def _SchoenhardtSurface(twist=numpy.pi / 6.0):
    """Schoenhardt polyhedron: the twisted triangular antiprism that admits NO
    tetrahedralization from its own six vertices - the classical solid whose
    boundary recovery needs an interior Steiner point."""
    angles = 2.0 * numpy.pi * numpy.arange(3) / 3.0
    bottom = numpy.c_[numpy.cos(angles), numpy.sin(angles), numpy.zeros(3)]
    top = numpy.c_[numpy.cos(angles + twist), numpy.sin(angles + twist), numpy.ones(3)]
    points = numpy.ascontiguousarray(numpy.vstack([bottom, top]))
    triangles = [[0, 2, 1], [3, 4, 5]]                     # caps, outward-wound
    for k in range(3):
        k1 = (k + 1) % 3
        triangles += [[k, k1, 3 + k1], [k, 3 + k1, 3 + k]]  # the reflex side split
    return points, numpy.array(triangles)


def _AsMesh(points, cells):
    return physicsnemo.mesh.Mesh(points=torch.tensor(points), cells=torch.tensor(cells))


@KratosUnittest.skipUnless(have_generate, _MISSING)
class TestFillSurfaceWithTetrahedra(KratosUnittest.TestCase):
    """The 3D interior fill upstream's fill_interior does not implement yet."""

    def _Fill(self, points, cells, extra=""):
        settings = Kratos.Parameters('{"full_output": true%s}' % extra)
        return generate.FillSurfaceWithTetrahedra(_AsMesh(points, cells), settings)

    def test_ConvexSolidIsFilledExactly(self):
        points, triangles = _CubeSurface()
        mesh, diagnostics = self._Fill(points, triangles)
        self.assertEqual(int(mesh.cells.shape[1]), 4)
        self.assertAlmostEqual(_TetVolume(mesh), 1.0, places=9)
        self.assertAlmostEqual(diagnostics["volume_ratio"], 1.0, places=9)
        self.assertAlmostEqual(diagnostics["boundary_area_ratio"], 1.0, places=9)
        self.assertTrue(diagnostics["boundary_manifold"])

    def test_NonConvexSolidIsFilledExactly(self):
        # The case a convex-hull tetrahedralization gets wrong without the
        # winding-number carve: the L's reentrant corner.
        points, triangles = _LPrismSurface()
        mesh, diagnostics = self._Fill(points, triangles)
        self.assertAlmostEqual(_TetVolume(mesh), 3.0, places=9)
        self.assertAlmostEqual(diagnostics["volume_ratio"], 1.0, places=9)
        self.assertAlmostEqual(diagnostics["boundary_area_ratio"], 1.0, places=9)
        self.assertTrue(diagnostics["boundary_manifold"])

    def test_InputVerticesSurviveBitIdentically(self):
        points, triangles = _CubeSurface()
        mesh, _ = self._Fill(points, triangles)
        self.assertTrue(torch.equal(
            torch.as_tensor(mesh.points)[:len(points)].double(),
            torch.tensor(points).double()))

    def test_EveryTetrahedronIsPositivelyOriented(self):
        # Delaunay does not wind consistently; Kratos needs positive Jacobians.
        for points, triangles in (_CubeSurface(), _LPrismSurface()):
            mesh, _ = self._Fill(points, triangles)
            vertices = torch.as_tensor(mesh.points).double()
            cells = torch.as_tensor(mesh.cells).long()
            a, b, c, d = (vertices[cells[:, i]] for i in range(4))
            signed = torch.einsum("ij,ij->i", b - a, torch.cross(c - a, d - a, dim=1))
            self.assertGreaterEqual(float(signed.min()), 0.0)

    def test_CurvedSurfaceMatchesItsOwnEnclosedVolume(self):
        n = 32
        axis = torch.linspace(-1.0, 1.0, n, dtype=torch.float64)
        x, y, z = torch.meshgrid(axis, axis, axis, indexing="ij")
        surface = generate.SurfaceFromLevelSet((x * x + y * y + z * z).sqrt() - 0.6)
        mesh, diagnostics = generate.FillSurfaceWithTetrahedra(
            surface, Kratos.Parameters('{"full_output": true}'))
        # The oracle is the polyhedron's own volume, not 4/3 pi r^3: the input
        # is a marching-cubes approximation of the sphere. The tolerance is
        # relative because retriangulating a non-planar boundary patch along
        # the other diagonal genuinely moves the enclosed volume a little
        # (~2e-5 here) - exactness is only available for planar-faced input.
        enclosed = abs(_SignedVolume(surface))
        self.assertAlmostEqual(_TetVolume(mesh), enclosed, delta=1e-3 * enclosed)
        self.assertTrue(diagnostics["boundary_manifold"])
        self.assertAlmostEqual(diagnostics["boundary_area_ratio"], 1.0, places=4)

    def test_UpstreamMethodStillReportsTheGap(self):
        points, triangles = _CubeSurface()
        with self.assertRaisesRegex(ValueError, "does not implement"):
            generate.FillSurfaceWithTetrahedra(
                _AsMesh(points, triangles), Kratos.Parameters('{"method": "upstream"}'))

    def test_OpenSurfaceIsRejectedUnlessStrictIsOff(self):
        points, triangles = _CubeSurface()
        open_surface = triangles[:-2]          # drop a face: no longer watertight
        with self.assertRaisesRegex(ValueError, "could not fill this surface"):
            generate.FillSurfaceWithTetrahedra(_AsMesh(points, open_surface))
        mesh, diagnostics = generate.FillSurfaceWithTetrahedra(
            _AsMesh(points, open_surface),
            Kratos.Parameters('{"strict": false, "full_output": true}'))
        self.assertEqual(int(mesh.cells.shape[1]), 4)
        self.assertIn("volume_ratio", diagnostics)

    def test_Validation(self):
        points, triangles = _CubeSurface()
        with self.assertRaisesRegex(ValueError, "method"):
            generate.FillSurfaceWithTetrahedra(
                _AsMesh(points, triangles), Kratos.Parameters('{"method": "netgen"}'))
        loop = physicsnemo.mesh.Mesh(
            points=torch.tensor([[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]]),
            cells=torch.tensor([[0, 1], [1, 2], [2, 0]]))
        with self.assertRaisesRegex(ValueError, "TRIANGLES"):
            generate.FillSurfaceWithTetrahedra(loop)

    def test_SchoenhardtSolidDefeatsTheDelaunayRoute(self):
        # The documented limitation, pinned: a solid needing Steiner points
        # for boundary recovery cannot be carved from its own points.
        points, triangles = _SchoenhardtSurface()
        with self.assertRaisesRegex(ValueError, "could not fill this surface"):
            generate.FillSurfaceWithTetrahedra(
                _AsMesh(points, triangles), Kratos.Parameters('{"method": "delaunay"}'))

    def test_FillModelPartWithTetrahedraIsSolvableKratos(self):
        points, triangles = _CubeSurface()
        model = Kratos.Model()
        skin = model.CreateModelPart("Skin")
        properties = skin.CreateNewProperties(1)
        for index, point in enumerate(points):
            skin.CreateNewNode(index + 1, float(point[0]), float(point[1]), float(point[2]))
        for index, triangle in enumerate(triangles):
            skin.CreateNewCondition("SurfaceCondition3D3N", index + 1,
                                    [int(node) + 1 for node in triangle], properties)

        volume = generate.FillModelPartWithTetrahedra(model, "Volume", skin)
        self.assertEqual(volume.NumberOfNodes(), len(points))
        self.assertGreater(volume.NumberOfElements(), 0)
        self.assertEqual(volume.ProcessInfo[Kratos.DOMAIN_SIZE], 3)
        self.assertEqual(min(node.Id for node in volume.Nodes), 1)
        for element in volume.Elements:
            self.assertEqual(len(element.GetGeometry()), 4)


@KratosUnittest.skipUnless(have_generate and have_tetgen,
                           "Missing required python modules: torch, physicsnemo, tetgen.")
class TestTetgenFill(KratosUnittest.TestCase):
    """The opt-in "tetgen" method: exact boundary recovery - input facets kept
    verbatim, Schoenhardt-class solids filled via interior Steiner points -
    which the Delaunay route documents as its two limitations."""

    def _Fill(self, points, cells, extra=""):
        settings = Kratos.Parameters('{"method": "tetgen", "full_output": true%s}' % extra)
        return generate.FillSurfaceWithTetrahedra(_AsMesh(points, cells), settings)

    def test_CubeKeepsAllTwelveFacets(self):
        # The Delaunay route's own docstring example: it keeps 8 of 12.
        points, triangles = _CubeSurface()
        mesh, diagnostics = self._Fill(points, triangles)
        self.assertTrue(diagnostics["facets_preserved"])
        self.assertAlmostEqual(_TetVolume(mesh), 1.0, places=9)
        self.assertAlmostEqual(diagnostics["volume_ratio"], 1.0, places=9)
        self.assertAlmostEqual(diagnostics["boundary_area_ratio"], 1.0, places=9)
        self.assertTrue(diagnostics["boundary_manifold"])

    def test_SchoenhardtSolidIsFilled(self):
        points, triangles = _SchoenhardtSurface()
        mesh, diagnostics = self._Fill(points, triangles)
        enclosed = abs(_SignedVolume(_AsMesh(points, triangles)))
        self.assertAlmostEqual(_TetVolume(mesh), enclosed, places=9)
        self.assertTrue(diagnostics["facets_preserved"])
        self.assertTrue(diagnostics["boundary_manifold"])
        # the classical solid is untetrahedralizable from its own vertices,
        # so at least one interior Steiner point is what made this possible
        self.assertGreaterEqual(diagnostics["steiner_points"], 1)

    def test_InputVerticesSurviveWithSteinerAppended(self):
        points, triangles = _CubeSurface()
        mesh, diagnostics = self._Fill(points, triangles)
        self.assertTrue(torch.equal(
            torch.as_tensor(mesh.points)[:len(points)].double(),
            torch.tensor(points).double()))
        self.assertEqual(
            int(mesh.points.shape[0]), len(points) + diagnostics["steiner_points"])

    def test_EveryTetrahedronIsPositivelyOriented(self):
        for points, triangles in (_CubeSurface(), _SchoenhardtSurface()):
            mesh, _ = self._Fill(points, triangles)
            vertices = torch.as_tensor(mesh.points).double()
            cells = torch.as_tensor(mesh.cells).long()
            a, b, c, d = (vertices[cells[:, i]] for i in range(4))
            signed = torch.einsum("ij,ij->i", b - a, torch.cross(c - a, d - a, dim=1))
            self.assertGreater(float(signed.min()), 0.0)


@KratosUnittest.skipUnless(have_generate, _MISSING)
class TestTetgenMissingPackage(KratosUnittest.TestCase):
    def test_MissingTetgenPackageIsActionable(self):
        # Simulate the package's absence regardless of the environment, so the
        # actionable error is pinned even where tetgen is installed.
        import sys

        class _Blocker:
            def find_spec(self, name, path=None, target=None):
                if name.split(".")[0] == "tetgen":
                    raise ImportError(f"'{name}' is blocked for this test")

        saved = {name: sys.modules.pop(name) for name in list(sys.modules)
                 if name.split(".")[0] == "tetgen"}
        sys.meta_path.insert(0, blocker := _Blocker())
        try:
            points, triangles = _CubeSurface()
            with self.assertRaisesRegex(ImportError, "pip install tetgen"):
                generate.FillSurfaceWithTetrahedra(
                    _AsMesh(points, triangles), Kratos.Parameters('{"method": "tetgen"}'))
        finally:
            sys.meta_path.remove(blocker)
            sys.modules.update(saved)


if __name__ == '__main__':
    KratosUnittest.main()
