"""Tests for surrogate-error-driven adaptive remeshing (size fields, MMG
volume adaptation, error-weighted surface partitioning)."""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import adaptive_remeshing
from KratosMultiphysics.PhysicsNeMoApplication.processes import adaptive_remesh_process
sys.path.insert(0, str(Path(__file__).parent / "kratos_solver_cases"))

have_meshing = kratos_utils.CheckIfApplicationsAvailable("MeshingApplication")
have_mmg = False
if have_meshing:
    import KratosMultiphysics.MeshingApplication as MeshingApplication
    have_mmg = hasattr(MeshingApplication, "MmgProcess2D")

have_convection_diffusion = kratos_utils.CheckIfApplicationsAvailable(
    "ConvectionDiffusionApplication", "LinearSolversApplication")

try:
    import torch
    import physicsnemo.mesh
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False

# remesh is Warp-backed since physicsnemo 2.2 (pyacvd is no longer involved)
have_remesh = have_physicsnemo


def _CreateTriangleSquare(model, name, divisions):
    """Unit-square Element2D3N mesh with `divisions` cells per axis."""
    model_part = model.CreateModelPart(name)
    model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
    model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
    props = model_part.CreateNewProperties(1)
    n = divisions + 1
    for i in range(n):
        for j in range(n):
            model_part.CreateNewNode(i * n + j + 1, i / divisions, j / divisions, 0.0)
    element_id = 0
    for i in range(divisions):
        for j in range(divisions):
            a, b = i * n + j + 1, (i + 1) * n + j + 1
            c, d = (i + 1) * n + j + 2, i * n + j + 2
            for triangle in ((a, b, c), (a, c, d)):
                element_id += 1
                model_part.CreateNewElement("Element2D3N", element_id, list(triangle), props)
    return model_part


class TestTargetSizeField(KratosUnittest.TestCase):
    """Pure Kratos: no torch, physicsnemo or MeshingApplication needed."""

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = _CreateTriangleSquare(self.model, "SizeField", 8)

    def test_EquidistributionMonotoneAndClipped(self):
        n = self.model_part.NumberOfNodes()
        error = numpy.full(n, 1e-3)
        error[0] = 1e-1   # far above target -> refine
        error[1] = 1e-9   # far below target -> coarsen (clipped)
        sizes = adaptive_remeshing.ComputeTargetSizeField(
            self.model_part, error, Kratos.Parameters("""{
                "target_error": 1e-3, "exponent": 0.5,
                "minimal_size": 1e-3, "maximal_size": 0.5
            }"""))
        self.assertEqual(sizes.shape, (n,))
        baseline = sizes[2]
        self.assertLess(sizes[0], baseline)       # high error refines
        self.assertGreater(sizes[1], baseline)    # low error coarsens
        self.assertGreaterEqual(sizes.min(), 1e-3)
        self.assertLessEqual(sizes.max(), 0.5)
        # at exactly the target error the current size is kept
        current_h = numpy.array([node.GetValue(Kratos.NODAL_H) for node in self.model_part.Nodes])
        numpy.testing.assert_allclose(sizes[2:], numpy.clip(current_h, 1e-3, 0.5)[2:], rtol=1e-12)

    def test_NodalErrorArrayCollapsesPerNode(self):
        residuals = {(1, "TEMPERATURE"): 0.5, (1, "VELOCITY_X"): 2.0, (2, "TEMPERATURE"): 0.1}
        collapsed = adaptive_remeshing.NodalErrorArray(self.model_part, residuals)
        self.assertEqual(collapsed.shape, (self.model_part.NumberOfNodes(),))
        self.assertAlmostEqual(collapsed[0], 2.0)   # node 1: max over DOFs
        self.assertAlmostEqual(collapsed[1], 0.1)
        self.assertAlmostEqual(collapsed[2], 0.0)   # unlisted node

    def test_ShapeMismatchRejected(self):
        with self.assertRaisesRegex(ValueError, "nodes"):
            adaptive_remeshing.ComputeTargetSizeField(
                self.model_part, numpy.ones(3), Kratos.Parameters("{}"))


@KratosUnittest.skipUnless(have_mmg, "MeshingApplication with MMG support is not available.")
class TestMmgAdaptation(KratosUnittest.TestCase):

    def test_LocalizedSizeFieldRefinesLocally(self):
        model = Kratos.Model()
        model_part = _CreateTriangleSquare(model, "Adapt", 10)
        for node in model_part.Nodes:
            node.SetSolutionStepValue(Kratos.TEMPERATURE, node.X + node.Y)
        nodes_before = model_part.NumberOfNodes()

        sizes = numpy.array([0.02 if node.X < 0.3 else 0.2 for node in model_part.Nodes])
        adaptive_remeshing.RunMmgAdaptation(model_part, sizes)

        self.assertTrue(model_part.Is(Kratos.MODIFIED))
        self.assertGreater(model_part.NumberOfNodes(), nodes_before)
        refined = [node.GetValue(Kratos.NODAL_H) for node in model_part.Nodes if node.X < 0.25]
        coarse = [node.GetValue(Kratos.NODAL_H) for node in model_part.Nodes if node.X > 0.5]
        self.assertLess(numpy.mean(refined), 0.5 * numpy.mean(coarse))
        # nodal values were interpolated onto the new mesh
        for node in model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE), node.X + node.Y, places=6)

    def test_SizeFieldShapeValidated(self):
        model = Kratos.Model()
        model_part = _CreateTriangleSquare(model, "AdaptBad", 4)
        with self.assertRaisesRegex(ValueError, "nodes"):
            adaptive_remeshing.RunMmgAdaptation(model_part, numpy.ones(3))


@KratosUnittest.skipUnless(have_mmg and have_convection_diffusion,
                           "Missing MeshingApplication (MMG) or ConvectionDiffusionApplication.")
class TestAdaptiveRemeshProcess(KratosUnittest.TestCase):

    def test_ResidualDrivenRemeshOnRealCase(self):
        import thermal_case

        model = Kratos.Model()
        analysis = thermal_case.CreateThermalAnalysis(model, conductivity=1.0, heat_flux=1.0)
        analysis.Run()
        model_part = model["ThermalModelPart"]

        # perturb the solved field so the residual concentrates in a corner
        for node in model_part.Nodes:
            if node.X < 0.3 and node.Y < 0.3:
                node.SetSolutionStepValue(
                    Kratos.TEMPERATURE, node.GetSolutionStepValue(Kratos.TEMPERATURE) + 0.5)

        nodes_before = model_part.NumberOfNodes()
        process = adaptive_remesh_process.Factory(Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "ThermalModelPart",
                "remesh_interval" : 1,
                "size_settings"   : { "target_error": 1e-4, "exponent": 0.5,
                                      "minimal_size": 0.02, "maximal_size": 0.2 }
            }
        }"""), model)
        model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        self.assertNotEqual(model_part.NumberOfNodes(), nodes_before)
        self.assertTrue(model_part.Is(Kratos.MODIFIED))

    def test_IntervalGating(self):
        import thermal_case

        model = Kratos.Model()
        analysis = thermal_case.CreateThermalAnalysis(model, conductivity=1.0, heat_flux=1.0)
        analysis.Run()
        model_part = model["ThermalModelPart"]
        nodes_before = model_part.NumberOfNodes()

        process = adaptive_remesh_process.Factory(Kratos.Parameters("""{
            "Parameters": { "model_part_name": "ThermalModelPart", "remesh_interval": 5 }
        }"""), model)
        model_part.ProcessInfo[Kratos.STEP] = 3  # not a multiple of 5
        process.ExecuteFinalizeSolutionStep()
        self.assertEqual(model_part.NumberOfNodes(), nodes_before)


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python modules: torch, physicsnemo.")
class TestSurfacePartition(KratosUnittest.TestCase):

    def _PlaneMesh(self, n=9):
        axis = torch.linspace(0.0, 1.0, n, dtype=torch.float32)
        xy = torch.stack(torch.meshgrid(axis, axis, indexing="ij"), dim=-1).reshape(-1, 2)
        points = torch.cat([xy, torch.zeros(len(xy), 1)], dim=1)
        cells = []
        for i in range(n - 1):
            for j in range(n - 1):
                a, b = i * n + j, (i + 1) * n + j
                c, d = (i + 1) * n + j + 1, i * n + j + 1
                cells += [[a, b, c], [a, c, d]]
        return physicsnemo.mesh.Mesh(points=points, cells=torch.tensor(cells, dtype=torch.int64))

    def test_WeightedSeedsFollowTheError(self):
        mesh = self._PlaneMesh()
        centroids = mesh.points[mesh.cells].mean(dim=1)
        weights = torch.where(centroids[:, 0] < 0.5, 50.0, 1.0)

        partition, seeds = adaptive_remeshing.WeightedSurfacePartition(
            mesh, 24, weights=weights, seed=0)
        left_fraction = float((seeds[:, 0] < 0.5).float().mean())
        self.assertGreater(left_fraction, 0.7)

        # every cell assigned, total area conserved
        self.assertEqual(tuple(partition.assignments.shape), (int(mesh.cells.shape[0]),))
        self.assertTrue(bool((partition.assignments >= 0).all()))
        self.assertAlmostEqual(float(partition.cluster_areas.sum()), 1.0, places=5)

    def test_UniformPartitionAndValidation(self):
        mesh = self._PlaneMesh(5)
        partition, seeds = adaptive_remeshing.WeightedSurfacePartition(mesh, 6, seed=0)
        self.assertEqual(tuple(seeds.shape), (6, 3))
        self.assertEqual(int(partition.assignments.max()) + 1 <= 6, True)
        with self.assertRaisesRegex(ValueError, "n_clusters"):
            adaptive_remeshing.WeightedSurfacePartition(mesh, 0)
        with self.assertRaisesRegex(ValueError, "weights"):
            adaptive_remeshing.WeightedSurfacePartition(mesh, 4, weights=torch.ones(3))


@KratosUnittest.skipUnless(have_remesh, "Missing required python modules: torch, physicsnemo.")
class TestRemeshSurface(KratosUnittest.TestCase):

    def test_IsotropicSurfaceRemesh(self):
        from physicsnemo.mesh.primitives.surfaces import sphere_icosahedral
        mesh = sphere_icosahedral.load(subdivisions=3)
        remeshed = adaptive_remeshing.RemeshSurface(mesh, 100)
        # since physicsnemo 2.2 the count targets output VERTICES exactly
        # (it targeted cells before), so this is an equality, not a band
        self.assertEqual(int(remeshed.points.shape[0]), 100)
        self.assertLess(int(remeshed.cells.shape[0]), int(mesh.cells.shape[0]))
        self.assertEqual(int(remeshed.cells.shape[1]), 3)
        # still a sphere: every vertex near unit radius
        radii = remeshed.points.norm(dim=1)
        self.assertLess(float((radii - 1.0).abs().max()), 0.1)


if __name__ == '__main__':
    KratosUnittest.main()
