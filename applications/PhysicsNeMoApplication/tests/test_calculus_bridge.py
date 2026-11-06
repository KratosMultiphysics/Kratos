"""Tests for the discrete-calculus bridge (physicsnemo.mesh.calculus on
tessellated Kratos meshes)."""

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges import calculus_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    import physicsnemo.mesh
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


def _MakePlaneSurfaceMesh(n=6):
    """Flat triangulated unit square at z=0, as a physicsnemo surface mesh."""
    axis = torch.linspace(0.0, 1.0, n, dtype=torch.float64)
    xy = torch.stack(torch.meshgrid(axis, axis, indexing="ij"), dim=-1).reshape(-1, 2)
    points = torch.cat([xy, torch.zeros(len(xy), 1, dtype=torch.float64)], dim=1)
    cells = []
    for i in range(n - 1):
        for j in range(n - 1):
            a, b = i * n + j, (i + 1) * n + j
            c, d = (i + 1) * n + j + 1, i * n + j + 1
            cells += [[a, b, c], [a, c, d]]
    return physicsnemo.mesh.Mesh(points=points, cells=torch.tensor(cells, dtype=torch.int64))


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python modules: torch, physicsnemo.")
class TestCalculusBridge(KratosUnittest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Calculus", 4,
            historical_variables=(Kratos.PRESSURE, Kratos.VELOCITY))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 2.0 * node.X + 3.0 * node.Y + 4.0 * node.Z)
            node.SetSolutionStepValue(Kratos.VELOCITY, [-node.Y, node.X, 0.0])

    def _BuildMesh(self, *specs):
        return domain_mesh_builder.BuildMesh(self.model_part, specs)

    def test_LsqGradientExactOnLinearField(self):
        mesh, _ = self._BuildMesh((Kratos.PRESSURE, "node_historical"))
        gradient = calculus_bridge.ComputeGradient(mesh, mesh.point_data["PRESSURE"])
        self.assertEqual(tuple(gradient.shape), (self.model_part.NumberOfNodes(), 3))
        expected = torch.tensor([2.0, 3.0, 4.0], dtype=gradient.dtype)
        self.assertLess(float((gradient - expected).abs().max()), 1e-8)

    def test_MultiChannelGradientIsChannelMajor(self):
        mesh, _ = self._BuildMesh()
        coordinates = torch.as_tensor(mesh.points)
        # An ASYMMETRIC Jacobian is essential here: u = (2y, 3z, 5x) gives
        # du_c/dx_d = [[0,2,0],[0,0,3],[5,0,0]], whose transpose differs. A
        # diagonal field (the obvious choice) is symmetric, so it cannot tell
        # the bridge's channel-major contract from upstream's derivative-first
        # layout - and upstream flipped that layout in physicsnemo 2.2.
        values = torch.stack([2.0 * coordinates[:, 1],
                              3.0 * coordinates[:, 2],
                              5.0 * coordinates[:, 0]], dim=1)
        gradient = calculus_bridge.ComputeGradient(mesh, values)
        self.assertEqual(tuple(gradient.shape), (self.model_part.NumberOfNodes(), 3, 3))
        expected = torch.tensor([[0.0, 2.0, 0.0],
                                 [0.0, 0.0, 3.0],
                                 [5.0, 0.0, 0.0]], dtype=gradient.dtype)
        self.assertLess(float((gradient - expected).abs().max()), 1e-8)

    def test_GradientLayoutIsIndependentOfChannelCount(self):
        """A non-square (C != D) field pins the axes unambiguously."""
        mesh, _ = self._BuildMesh()
        coordinates = torch.as_tensor(mesh.points)
        values = torch.stack([coordinates[:, 0], 4.0 * coordinates[:, 1]], dim=1)  # C=2, D=3
        gradient = calculus_bridge.ComputeGradient(mesh, values)
        self.assertEqual(tuple(gradient.shape), (self.model_part.NumberOfNodes(), 2, 3))
        expected = torch.tensor([[1.0, 0.0, 0.0], [0.0, 4.0, 0.0]], dtype=gradient.dtype)
        self.assertLess(float((gradient - expected).abs().max()), 1e-8)

    def test_ScalarGradientKeepsItsShape(self):
        mesh, _ = self._BuildMesh((Kratos.PRESSURE, "node_historical"))
        gradient = calculus_bridge.ComputeGradient(mesh, mesh.point_data["PRESSURE"])
        self.assertEqual(tuple(gradient.shape), (self.model_part.NumberOfNodes(), 3))

    def test_DecGradientRefusedOnVolumeMeshes(self):
        mesh, _ = self._BuildMesh((Kratos.PRESSURE, "node_historical"),
                                  (Kratos.VELOCITY, "node_historical"))
        with self.assertRaisesRegex(ValueError, "volume"):
            calculus_bridge.ComputeGradient(mesh, mesh.point_data["PRESSURE"], method="dec")
        with self.assertRaisesRegex(ValueError, "volume"):
            calculus_bridge.ComputeDivergence(mesh, mesh.point_data["VELOCITY"], method="dec")

    def test_CurlAndDivergenceOfRigidRotation(self):
        mesh, _ = self._BuildMesh((Kratos.VELOCITY, "node_historical"))
        velocity = mesh.point_data["VELOCITY"]
        curl = calculus_bridge.ComputeCurl(mesh, velocity)
        expected = torch.tensor([0.0, 0.0, 2.0], dtype=curl.dtype)
        self.assertLess(float((curl - expected).abs().max()), 1e-8)
        divergence = calculus_bridge.ComputeDivergence(mesh, velocity)
        self.assertLess(float(divergence.abs().max()), 1e-8)

    def test_LaplacianInteriorZeroOnLinearField(self):
        mesh, _ = self._BuildMesh((Kratos.PRESSURE, "node_historical"))
        laplacian = calculus_bridge.ComputeLaplacian(mesh, mesh.point_data["PRESSURE"])
        mask = calculus_bridge.InteriorPointMask(mesh)
        self.assertGreater(int(mask.sum()), 0)
        self.assertLess(int(mask.sum()), self.model_part.NumberOfNodes())
        self.assertLess(float(laplacian[mask].abs().max()), 1e-6)

    def test_InteriorPointMaskMatchesLattice(self):
        mesh, _ = self._BuildMesh()
        mask = calculus_bridge.InteriorPointMask(mesh).numpy()
        # structured unit-cube lattice: interior nodes have no 0/1 coordinate
        coordinates = numpy.asarray(mesh.points)
        on_boundary = numpy.any((coordinates <= 1e-12) | (coordinates >= 1.0 - 1e-12), axis=1)
        numpy.testing.assert_array_equal(mask, ~on_boundary)

    def test_IntegrateFieldGivesVolume(self):
        mesh, _ = self._BuildMesh()
        ones = torch.ones(self.model_part.NumberOfNodes(), dtype=torch.float64)
        volume = calculus_bridge.IntegrateField(mesh, ones, data_source="points")
        self.assertAlmostEqual(float(volume), 1.0, places=10)

    def test_SurfaceIntrinsicGradientAndChannelLoop(self):
        mesh = _MakePlaneSurfaceMesh()
        x, y = mesh.points[:, 0], mesh.points[:, 1]
        scalar = 2.0 * x + 3.0 * y
        gradient = calculus_bridge.ComputeGradient(mesh, scalar)
        expected = torch.tensor([2.0, 3.0, 0.0], dtype=gradient.dtype)
        self.assertLess(float((gradient - expected).abs().max()), 1e-8)
        # multi-channel intrinsic path (crashed upstream before 2.2)
        stacked = torch.stack([scalar, -scalar], dim=1)
        gradient2 = calculus_bridge.ComputeGradient(mesh, stacked)
        self.assertEqual(tuple(gradient2.shape), (mesh.points.shape[0], 2, 3))
        self.assertLess(float((gradient2[:, 0] - expected).abs().max()), 1e-8)
        self.assertLess(float((gradient2[:, 1] + expected).abs().max()), 1e-8)

    def test_GradientIsDifferentiable(self):
        mesh, _ = self._BuildMesh()
        values = torch.as_tensor(mesh.points)[:, 0].clone().requires_grad_(True)
        gradient = calculus_bridge.ComputeGradient(mesh, values)
        gradient.square().sum().backward()
        self.assertIsNotNone(values.grad)
        self.assertGreater(float(values.grad.abs().sum()), 0.0)

    def test_ComputeNodalDerivativesWritesBack(self):
        written = calculus_bridge.ComputeNodalDerivatives(self.model_part, Kratos.Parameters("""{
            "operations": [
                { "field": "PRESSURE", "operation": "gradient",
                  "output_variable": "PRESSURE_GRADIENT" },
                { "field": "VELOCITY", "operation": "divergence",
                  "output_variable": "NODAL_ERROR" }
            ]
        }"""))
        self.assertEqual(set(written), {"PRESSURE_GRADIENT", "NODAL_ERROR"})
        for node in self.model_part.Nodes:
            gradient = node.GetValue(Kratos.PRESSURE_GRADIENT)
            self.assertAlmostEqual(gradient[0], 2.0, places=8)
            self.assertAlmostEqual(gradient[1], 3.0, places=8)
            self.assertAlmostEqual(gradient[2], 4.0, places=8)
            self.assertAlmostEqual(node.GetValue(Kratos.NODAL_ERROR), 0.0, places=8)

    def test_ComputeNodalDerivativesZeroBoundary(self):
        calculus_bridge.ComputeNodalDerivatives(self.model_part, Kratos.Parameters("""{
            "operations": [
                { "field": "PRESSURE", "operation": "laplacian",
                  "output_variable": "NODAL_ERROR" }
            ],
            "zero_boundary": true
        }"""))
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(node.GetValue(Kratos.NODAL_ERROR), 0.0, places=5)

    def test_ComputeNodalDerivativesValidation(self):
        with self.assertRaisesRegex(ValueError, "Unknown operation"):
            calculus_bridge.ComputeNodalDerivatives(self.model_part, Kratos.Parameters("""{
                "operations": [ { "field": "PRESSURE", "operation": "hessian",
                                  "output_variable": "NODAL_ERROR" } ]
            }"""))

    def test_IntegrateNodalField(self):
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0)
        self.assertAlmostEqual(
            calculus_bridge.IntegrateNodalField(self.model_part, Kratos.PRESSURE), 1.0,
            places=10)


if __name__ == '__main__':
    KratosUnittest.main()
