"""Tests for the differentiable shape-deformation layer (physicsnemo >= 2.2
deformers and energies) and the shape sensitivities built on it."""

import numpy

import sys
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
from KratosMultiphysics.PhysicsNeMoApplication.physics import differentiable_residual
from KratosMultiphysics.PhysicsNeMoApplication.physics import sensitivity_utils
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import deformation
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    from physicsnemo.mesh import deformation as _pn_deformation  # noqa: F401
    have_deformation = True
except ImportError:
    have_deformation = False

_MISSING = "Missing required python modules: torch, physicsnemo >= 2.2."

have_mesh_moving = kratos_utils.CheckIfApplicationsAvailable(
    "MeshMovingApplication", "LinearSolversApplication")
have_convection_diffusion = kratos_utils.CheckIfApplicationsAvailable(
    "ConvectionDiffusionApplication", "LinearSolversApplication")

sys.path.insert(0, str(Path(__file__).parent / "kratos_solver_cases"))


@KratosUnittest.skipUnless(have_deformation, _MISSING)
class TestDeformPoints(KratosUnittest.TestCase):

    def setUp(self):
        torch.manual_seed(0)
        self.points = torch.rand(40, 3, dtype=torch.float64)

    def test_FfdUniformControlIsExactTranslation(self):
        shift = torch.tensor([0.1, 0.2, -0.3], dtype=torch.float64)
        lattice = shift.expand(2, 2, 2, 3).clone()
        deformed = deformation.DeformPoints(
            self.points, lattice, "ffd", origin=[0.0, 0.0, 0.0], extent=[1.0, 1.0, 1.0])
        self.assertLess(float((deformed - (self.points + shift)).abs().max()), 1e-10)

    def test_ZeroControlIsTheIdentity(self):
        for method, controls, options in (
                ("ffd", torch.zeros(2, 2, 2, 3, dtype=torch.float64), {}),
                ("displace", torch.zeros_like(self.points), {}),
        ):
            deformed = deformation.DeformPoints(self.points, controls, method, **options)
            self.assertLess(float((deformed - self.points).abs().max()), 1e-12,
                            msg=f"{method} is not the identity at zero control")

    def test_RbfInterpolatesItsControlDisplacements(self):
        torch.manual_seed(1)
        # control points must be in general position: a degenerate (e.g.
        # coplanar) set makes the thin-plate-spline system singular upstream
        control_points = torch.rand(8, 3, dtype=torch.float64) * 2.0 - 0.5
        control_displacements = torch.rand(8, 3, dtype=torch.float64) * 0.05
        moved = deformation.DeformPoints(
            control_points, control_displacements, "rbf", control_points=control_points)
        self.assertLess(
            float((moved - (control_points + control_displacements)).abs().max()), 1e-6)

    def test_DisplaceIsPerPoint(self):
        displacement = torch.zeros_like(self.points)
        displacement[0] = torch.tensor([1.0, 0.0, 0.0], dtype=torch.float64)
        deformed = deformation.DeformPoints(self.points, displacement, "displace")
        self.assertAlmostEqual(float(deformed[0, 0] - self.points[0, 0]), 1.0, places=12)
        self.assertLess(float((deformed[1:] - self.points[1:]).abs().max()), 1e-12)

    def test_Validation(self):
        with self.assertRaisesRegex(ValueError, "Unknown deformation method"):
            deformation.DeformPoints(self.points, torch.zeros(3), "spline")
        with self.assertRaisesRegex(ValueError, "control_points"):
            deformation.DeformPoints(self.points, torch.zeros(4, 3), "rbf")

    def test_GradientsFlowToControlParameters(self):
        controls = torch.zeros(2, 2, 2, 3, dtype=torch.float64, requires_grad=True)
        deformed = deformation.DeformPoints(
            self.points, controls, "ffd", origin=[0.0, 0.0, 0.0], extent=[1.0, 1.0, 1.0])
        deformed.sum().backward()
        self.assertIsNotNone(controls.grad)
        self.assertGreater(float(controls.grad.abs().sum()), 0.0)


@KratosUnittest.skipUnless(have_deformation, _MISSING)
class TestRegularizationEnergy(KratosUnittest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Energy", 3, historical_variables=(Kratos.PRESSURE,))
        self.mesh, _ = domain_mesh_builder.BuildMesh(self.model_part)
        self.points = torch.as_tensor(self.mesh.points)

    def test_StrainVanishesAtReferenceAndGrowsWithDistortion(self):
        at_rest = deformation.RegularizationEnergy(self.mesh, self.points, "strain")
        self.assertAlmostEqual(float(at_rest), 0.0, places=9)
        stretched = self.points * torch.tensor([1.3, 1.0, 1.0], dtype=self.points.dtype)
        distorted = deformation.RegularizationEnergy(self.mesh, stretched, "strain")
        self.assertGreater(float(distorted), 0.0)

    def test_InversionEnergyDetectsFoldedElements(self):
        healthy = deformation.RegularizationEnergy(self.mesh, self.points, "inversion")
        folded = self.points.clone()
        folded[:, 0] = -folded[:, 0]   # mirror: every tet's orientation flips
        inverted = deformation.RegularizationEnergy(self.mesh, folded, "inversion")
        self.assertGreater(float(inverted), float(healthy))

    def test_EnergyIsDifferentiableWrtPoints(self):
        points = self.points.clone().requires_grad_(True)
        stretched = points * 1.2
        deformation.RegularizationEnergy(self.mesh, stretched, "strain").backward()
        self.assertIsNotNone(points.grad)
        self.assertGreater(float(points.grad.abs().sum()), 0.0)

    def test_UnknownEnergyRejected(self):
        with self.assertRaisesRegex(ValueError, "Unknown energy"):
            deformation.RegularizationEnergy(self.mesh, self.points, "bogus")


@KratosUnittest.skipUnless(have_deformation, _MISSING)
class TestCoordinateWriteBack(KratosUnittest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Move", 2, historical_variables=(Kratos.DISPLACEMENT,))

    def test_WriteCoordinatesRoundTrips(self):
        original = graph_bridge.NodePositions(self.model_part)
        shifted = original + numpy.array([0.1, -0.2, 0.3])
        deformation.WriteNodeCoordinates(self.model_part, shifted)
        numpy.testing.assert_allclose(
            graph_bridge.NodePositions(self.model_part), shifted, atol=1e-12)
        # X0 is the reference configuration and must NOT move
        reference = numpy.array([[n.X0, n.Y0, n.Z0] for n in self.model_part.Nodes])
        numpy.testing.assert_allclose(reference, original, atol=1e-12)

    def test_UpdateDisplacementMatchesTheOffset(self):
        original = graph_bridge.NodePositions(self.model_part)
        shifted = original + numpy.array([0.05, 0.0, 0.0])
        deformation.WriteNodeCoordinates(self.model_part, shifted, update_displacement=True)
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.DISPLACEMENT_X), 0.05, places=12)

    def test_RowCountValidated(self):
        with self.assertRaisesRegex(ValueError, "nodes"):
            deformation.WriteNodeCoordinates(self.model_part, numpy.zeros((3, 3)))

    def test_DeformModelPartMovesTheMesh(self):
        original = graph_bridge.NodePositions(self.model_part)
        shift = torch.tensor([0.02, 0.0, -0.01], dtype=torch.float64)
        coordinates = deformation.DeformModelPart(
            self.model_part, shift.expand(2, 2, 2, 3).clone(), "ffd")
        numpy.testing.assert_allclose(
            coordinates, original + shift.numpy(), atol=1e-10)
        numpy.testing.assert_allclose(
            graph_bridge.NodePositions(self.model_part), coordinates, atol=1e-12)


@KratosUnittest.skipUnless(have_deformation, _MISSING)
class TestShapeSensitivities(KratosUnittest.TestCase):
    """dJ/d(control) through a deformation and a surrogate, against FD."""

    def setUp(self):
        torch.manual_seed(0)
        self.points = torch.rand(30, 3, dtype=torch.float64)
        self.features = torch.ones(30, 1, dtype=torch.float64)
        # a tiny geometry-sensitive surrogate: it reads the coordinates
        self.model = torch.nn.Sequential(
            torch.nn.Linear(4, 8), torch.nn.Tanh(), torch.nn.Linear(8, 1)).double()

    def _Objective(self, prediction):
        return prediction.square().mean()

    def test_ShapeGradientMatchesFiniteDifferences(self):
        controls = torch.zeros(2, 2, 2, 3, dtype=torch.float64)
        options = {"origin": [0.0, 0.0, 0.0], "extent": [1.0, 1.0, 1.0]}
        result = sensitivity_utils.ComputeShapeSensitivities(
            self.model, self.features, controls, self._Objective, self.points,
            method="ffd", **options)

        self.assertEqual(result["control_displacements"].shape, (2, 2, 2, 3))
        self.assertEqual(result["coordinates"].shape, (30, 3))

        def objective_at(control_array):
            deformed = deformation.DeformPoints(
                self.points, torch.as_tensor(control_array), "ffd", **options)
            # mirror the "generic" interface exactly: coordinates PREPENDED to
            # the features, with a leading batch axis
            inputs = torch.cat([deformed, self.features], dim=-1)[None]
            with torch.no_grad():
                return float(self._Objective(self.model(inputs)))

        step = 1e-6
        base = controls.numpy().copy()
        for index in [(0, 0, 0, 0), (1, 0, 1, 2), (1, 1, 1, 1)]:
            plus, minus = base.copy(), base.copy()
            plus[index] += step
            minus[index] -= step
            finite_difference = (objective_at(plus) - objective_at(minus)) / (2.0 * step)
            self.assertAlmostEqual(
                result["control_displacements"][index], finite_difference, places=6,
                msg=f"adjoint vs FD mismatch at control {index}")

    def test_ZeroSensitivityForAShapeBlindObjective(self):
        constant = torch.nn.Sequential(torch.nn.Linear(4, 1)).double()
        with torch.no_grad():
            constant[0].weight.zero_()
            constant[0].bias.fill_(1.0)
        result = sensitivity_utils.ComputeShapeSensitivities(
            constant, self.features, torch.zeros(2, 2, 2, 3, dtype=torch.float64),
            self._Objective, self.points, method="ffd",
            origin=[0.0, 0.0, 0.0], extent=[1.0, 1.0, 1.0])
        self.assertLess(numpy.abs(result["control_displacements"]).max(), 1e-12)


@KratosUnittest.skipUnless(have_deformation, _MISSING)
@KratosUnittest.skipUnless(have_mesh_moving,
                           "Missing MeshMoving/LinearSolvers applications.")
class TestMeshMovingInteriorSmoothing(KratosUnittest.TestCase):
    """The shipped deformation layer driving MeshMovingApplication.

    The deformers move a *boundary*; the interior has to follow without
    folding. MeshMovingApplication is the Kratos-side counterpart that does
    that, so the pairing is the real deliverable: deform the boundary with
    DeformPoints, hand the boundary displacement to the mesh solver, and let
    it smooth the interior.
    """

    _DIVISIONS = 4
    _SHEAR = 0.35   # large enough that leaving the interior put folds elements

    def setUp(self):
        torch.manual_seed(0)
        self.model = Kratos.Model()
        self.model_part = self._CreateBoxMesh(self.model, self._DIVISIONS)

    @staticmethod
    def _CreateBoxMesh(model, divisions):
        model_part = model.CreateModelPart("MeshMotion")
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        model_part.AddNodalSolutionStepVariable(Kratos.MESH_DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(Kratos.MESH_REACTION)
        properties = model_part.CreateNewProperties(1)

        n = divisions
        node_id = lambda i, j, k: i * (n + 1) * (n + 1) + j * (n + 1) + k + 1
        for i in range(n + 1):
            for j in range(n + 1):
                for k in range(n + 1):
                    model_part.CreateNewNode(node_id(i, j, k), i / n, j / n, k / n)
        # Freudenthal subdivision: consistent, no hanging faces
        tets = [(0, 1, 3, 7), (0, 1, 7, 5), (0, 5, 7, 4),
                (0, 3, 2, 7), (0, 6, 4, 7), (0, 2, 6, 7)]
        corner = lambda i, j, k, c: node_id(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1))
        element_id = 0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for tet in tets:
                        element_id += 1
                        model_part.CreateNewElement(
                            "Element3D4N", element_id,
                            [corner(i, j, k, c) for c in tet], properties)
        model_part.SetBufferSize(2)
        return model_part

    @staticmethod
    def _IsOnBoundary(node):
        return (min(node.X0, node.Y0, node.Z0) <= 1e-12
                or max(node.X0, node.Y0, node.Z0) >= 1.0 - 1e-12)

    @staticmethod
    def _ElementVolumes(model_part):
        return numpy.array([element.GetGeometry().Volume()
                            for element in model_part.Elements])

    def _BoundaryDisplacementFromDeformer(self):
        """A shear of the top face, produced by the shipped FFD deformer."""
        boundary = [node for node in self.model_part.Nodes if self._IsOnBoundary(node)]
        points = torch.tensor([[node.X0, node.Y0, node.Z0] for node in boundary],
                              dtype=torch.float64)
        # an FFD lattice whose top layer is shifted in +x: a pure shear
        lattice = torch.zeros(2, 2, 2, 3, dtype=torch.float64)
        lattice[:, :, 1, 0] = self._SHEAR
        deformed = deformation.DeformPoints(
            points, lattice, "ffd", origin=[0.0, 0.0, 0.0], extent=[1.0, 1.0, 1.0])
        return boundary, (deformed - points).numpy()

    def _RunMeshMoving(self):
        from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import (
            MeshMovingAnalysis)

        parameters = Kratos.Parameters("""{
            "problem_data"    : {
                "problem_name" : "mesh_motion", "parallel_type" : "OpenMP",
                "echo_level"   : 0, "start_time" : 0.0, "end_time" : 1.0
            },
            "solver_settings" : {
                "solver_type"             : "structural_similarity",
                "model_part_name"         : "MeshMotion",
                "domain_size"             : 3,
                "echo_level"              : 0,
                "calculate_mesh_velocity" : false,
                "model_import_settings"   : { "input_type" : "use_input_model_part" },
                "time_stepping"           : { "time_step" : 1.0 }
            },
            "processes"       : {}
        }""")
        MeshMovingAnalysis(self.model, parameters).Run()

    def test_DeformedBoundaryIsSmoothedIntoTheInterior(self):
        boundary, displacement = self._BoundaryDisplacementFromDeformer()
        self.assertGreater(numpy.abs(displacement).max(), 0.1)  # a real deformation

        for node, offset in zip(boundary, displacement):
            node.SetSolutionStepValue(Kratos.MESH_DISPLACEMENT, [float(v) for v in offset])
            for component in (Kratos.MESH_DISPLACEMENT_X, Kratos.MESH_DISPLACEMENT_Y,
                              Kratos.MESH_DISPLACEMENT_Z):
                node.Fix(component)

        self._RunMeshMoving()

        # the boundary the deformer prescribed is reproduced exactly
        for node, offset in zip(boundary, displacement):
            moved = node.GetSolutionStepValue(Kratos.MESH_DISPLACEMENT)
            for value, reference in zip(moved, offset):
                self.assertAlmostEqual(value, reference, places=12)

        # the interior followed rather than staying put
        interior_moved = sum(
            1 for node in self.model_part.Nodes if not self._IsOnBoundary(node)
            and numpy.linalg.norm(node.GetSolutionStepValue(Kratos.MESH_DISPLACEMENT)) > 1e-9)
        interior_total = sum(1 for node in self.model_part.Nodes
                             if not self._IsOnBoundary(node))
        self.assertEqual(interior_moved, interior_total)

        # and the mesh is still usable
        volumes = self._ElementVolumes(self.model_part)
        self.assertGreater(volumes.min(), 0.0)

    def test_SmoothingBeatsMovingTheBoundaryAlone(self):
        # The negative control: applying the same boundary displacement with
        # the interior left behind is what mesh smoothing exists to avoid.
        boundary, displacement = self._BoundaryDisplacementFromDeformer()
        for node, offset in zip(boundary, displacement):
            node.X = node.X0 + float(offset[0])
            node.Y = node.Y0 + float(offset[1])
            node.Z = node.Z0 + float(offset[2])
        unsmoothed = self._ElementVolumes(self.model_part).min()

        # reset and do it properly
        for node in self.model_part.Nodes:
            node.X, node.Y, node.Z = node.X0, node.Y0, node.Z0
        for node, offset in zip(boundary, displacement):
            node.SetSolutionStepValue(Kratos.MESH_DISPLACEMENT, [float(v) for v in offset])
            for component in (Kratos.MESH_DISPLACEMENT_X, Kratos.MESH_DISPLACEMENT_Y,
                              Kratos.MESH_DISPLACEMENT_Z):
                node.Fix(component)
        self._RunMeshMoving()
        smoothed = self._ElementVolumes(self.model_part).min()

        self.assertGreater(smoothed, unsmoothed)


@KratosUnittest.skipUnless(have_deformation, _MISSING)
@KratosUnittest.skipUnless(have_convection_diffusion,
                           "Missing ConvectionDiffusion/LinearSolvers applications.")
class TestExactControlSensitivities(KratosUnittest.TestCase):
    """The FEM-exact chain: dJ/dX from the adjoint, pushed onto controls.

    ComputeShapeSensitivities differentiates a SURROGATE's forward pass;
    this path takes the discretely exact nodal gradient from
    ComputeShapeSensitivityField and applies only the deformation's chain
    rule, so the control gradient is as exact as the FEM adjoint behind it.
    """

    # a lattice spanning the unit square, with a non-degenerate z extent -
    # the planar mesh's own bounding box would give extent z = 0
    _ORIGIN = [0.0, 0.0, -0.5]
    _EXTENT = [1.0, 1.0, 1.0]

    @staticmethod
    def _SolveThermal(control=None, reference=None):
        """Solves the thermal case, optionally on an FFD-deformed mesh.

        The case data (including which nodes are Dirichlet) is applied by
        CreateThermalAnalysis on the UNDEFORMED mesh, so deforming after
        construction keeps the boundary set fixed - which is what the
        adjoint assumes.
        """
        import thermal_case
        model = Kratos.Model()
        analysis = thermal_case.CreateThermalAnalysis(
            model, conductivity=2.0, heat_flux=1.0, divisions=6)
        analysis.Initialize()
        model_part = model["ThermalModelPart"]

        if control is not None:
            points = torch.as_tensor(reference, dtype=torch.float64)
            deformed = deformation.DeformPoints(
                points, torch.as_tensor(control, dtype=torch.float64), "ffd",
                origin=TestExactControlSensitivities._ORIGIN,
                extent=TestExactControlSensitivities._EXTENT).numpy()
            # write BOTH configurations: which one an element reads is
            # element-dependent, and this is meant to be a different mesh
            for node, position in zip(model_part.Nodes, deformed):
                node.X0, node.Y0, node.Z0 = (float(position[0]), float(position[1]),
                                             float(position[2]))
                node.X, node.Y, node.Z = node.X0, node.Y0, node.Z0

        analysis.RunSolutionLoop()
        return model, model_part

    @staticmethod
    def _Objective(model_part):
        return sum(node.GetSolutionStepValue(Kratos.TEMPERATURE)
                   for node in model_part.Nodes)

    def test_EveryMethodReturnsAControlShapedGradient(self):
        rng = numpy.random.default_rng(0)
        reference = rng.random((20, 3))
        field = rng.random((20, 3))
        control_points = rng.random((5, 3))

        cases = (
            ("ffd", numpy.zeros((3, 3, 3, 3)), {"origin": [0.0, 0.0, 0.0],
                                                "extent": [1.0, 1.0, 1.0]}),
            ("rbf", numpy.zeros((5, 3)), {"control_points": control_points}),
            ("morph", numpy.zeros((5, 3)), {"control_points": control_points,
                                            "radius": 0.5}),
            ("displace", numpy.zeros((20, 3)), {}),
        )
        for method, control, options in cases:
            with self.subTest(method=method):
                gradient = sensitivity_utils.ComputeControlSensitivities(
                    field, reference, control, method, **options)
                self.assertEqual(gradient.shape, control.shape)
                self.assertTrue(numpy.all(numpy.isfinite(gradient)))

    def test_MismatchedFieldIsRejected(self):
        rng = numpy.random.default_rng(1)
        with self.assertRaisesRegex(ValueError, "row for row"):
            sensitivity_utils.ComputeControlSensitivities(
                rng.random((7, 3)), rng.random((20, 3)), numpy.zeros((20, 3)), "displace")

    def test_ControlGradientMatchesFiniteDifferencesThroughRealSolves(self):
        # The end-to-end claim: adjoint field -> deformer VJP -> dJ/d(control),
        # against re-solving the FEM on meshes deformed by the control itself.
        model, model_part = self._SolveThermal()
        reference = numpy.array([[node.X0, node.Y0, node.Z0]
                                 for node in model_part.Nodes])

        assembler = differentiable_residual.TangentAssembler(model_part)
        dof_map = differentiable_residual.DofFieldMap(
            assembler, [("TEMPERATURE", "node_historical")])
        dJ_du = numpy.ones(dof_map.n_equations)   # J = sum of temperatures

        field = sensitivity_utils.ComputeShapeSensitivityField(
            assembler, dof_map, dJ_du, fd_step=1e-6)

        control = numpy.zeros((2, 2, 2, 3))
        gradient = sensitivity_utils.ComputeControlSensitivities(
            field, reference, control, "ffd",
            origin=self._ORIGIN, extent=self._EXTENT)
        self.assertEqual(gradient.shape, control.shape)

        # a few lattice entries that actually move the mesh in-plane
        entries = ((1, 0, 0, 0), (1, 1, 0, 1), (0, 1, 1, 0))
        step = 1e-5
        checked = 0
        for entry in entries:
            plus = control.copy()
            plus[entry] += step
            minus = control.copy()
            minus[entry] -= step
            _, part_plus = self._SolveThermal(plus, reference)
            _, part_minus = self._SolveThermal(minus, reference)
            finite_difference = (self._Objective(part_plus)
                                 - self._Objective(part_minus)) / (2.0 * step)

            with self.subTest(entry=entry):
                self.assertGreater(abs(finite_difference), 1e-9)
                self.assertAlmostEqual(gradient[entry], finite_difference,
                                       delta=2e-4 * abs(finite_difference))
            checked += 1
        self.assertEqual(checked, len(entries))


if __name__ == '__main__':
    KratosUnittest.main()
