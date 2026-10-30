"""Tests for sensitivity_utils: surrogate autograd sensitivities and exact
adjoint parameter sensitivities (ConvectionDiffusion-gated), validated
against finite differences through real solves."""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.physics import differentiable_residual
from KratosMultiphysics.PhysicsNeMoApplication.physics import sensitivity_utils
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

have_convection_diffusion = kratos_utils.CheckIfApplicationsAvailable(
    "ConvectionDiffusionApplication", "LinearSolversApplication")

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))


def _SolveThermal(conductivity, heat_flux=1.0, divisions=4):
    import thermal_case
    model = Kratos.Model()
    analysis = thermal_case.CreateThermalAnalysis(
        model, conductivity=conductivity, heat_flux=heat_flux, divisions=divisions)
    analysis.Run()
    return model, model["ThermalModelPart"]


def _TotalTemperature(model_part):
    return sum(node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in model_part.Nodes)


@KratosUnittest.skipUnless(have_torch and have_convection_diffusion,
                           "Missing torch or ConvectionDiffusion/LinearSolvers applications.")
class TestSensitivityUtils(KratosUnittest.TestCase):
    def _AdjointSensitivity(self, parameter_variable, theta0, conductivity, heat_flux):
        """dJ/dtheta of J = sum(T) via the adjoint path, for a nodal scalar
        parameter applied uniformly."""
        model, model_part = _SolveThermal(conductivity, heat_flux)
        assembler = differentiable_residual.TangentAssembler(model_part)
        dof_map = differentiable_residual.DofFieldMap(
            assembler, [("TEMPERATURE", "node_historical")])

        # J = sum over nodes of T = sum over equations of u (fixed T = 0)
        dJ_du = numpy.ones(dof_map.n_equations)

        def apply_theta(value):
            for node in model_part.Nodes:
                node.SetSolutionStepValue(parameter_variable, value)

        sensitivities = sensitivity_utils.ComputeParameterSensitivities(
            assembler, dof_map, dJ_du,
            {"theta": (apply_theta, theta0)}, fd_step=1e-6)
        return sensitivities["theta"], _TotalTemperature(model_part)

    def test_AdjointMatchesFullFdAndAnalytics_Conductivity(self):
        k0 = 2.0
        adjoint, J0 = self._AdjointSensitivity(Kratos.CONDUCTIVITY, k0, k0, 1.0)

        # THE referee: full finite differences through two REAL solves -
        # the adjoint matches to ~8 significant digits
        h = 1e-5
        _, part_plus = _SolveThermal(k0 + h, 1.0)
        _, part_minus = _SolveThermal(k0 - h, 1.0)
        full_fd = (_TotalTemperature(part_plus) - _TotalTemperature(part_minus)) / (2.0 * h)
        self.assertAlmostEqual(adjoint, full_fd, delta=1e-6 * abs(full_fd))

        # the continuous analytic -J/k is only a LOOSE sanity bound: the
        # stabilized element's tau makes T ~ f/k approximately linear only
        self.assertAlmostEqual(adjoint, -J0 / k0, delta=0.05 * abs(J0 / k0))

    def test_AdjointLoadSensitivityExact(self):
        f0 = 3.0
        adjoint, J0 = self._AdjointSensitivity(Kratos.HEAT_FLUX, f0, 1.5, f0)
        # linear in the load: dJ/df = J/f, a hard equality
        self.assertAlmostEqual(adjoint, J0 / f0, delta=1e-6 * abs(J0 / f0))

    def test_SurrogateSensitivitiesFiniteDifference(self):
        torch.manual_seed(0)
        model = torch.nn.Sequential(
            torch.nn.Linear(4, 8), torch.nn.Tanh(), torch.nn.Linear(8, 1)).double()
        rng = numpy.random.default_rng(0)
        coordinates = rng.random((12, 3))
        features = rng.random((12, 1))

        def objective(prediction):
            return prediction.square().mean()

        result = sensitivity_utils.ComputeSurrogateSensitivities(
            model, features, coordinates, objective,
            model_interface="generic", wrt=("coordinates", "features"))
        self.assertEqual(result["coordinates"].shape, (12, 3))
        self.assertEqual(result["features"].shape, (12, 1))

        # central FD spot checks on a few coordinate entries
        def J_of(coords):
            with torch.no_grad():
                x = torch.tensor(numpy.concatenate([coords, features], axis=1))
                return float(model(x).square().mean())

        eps = 1e-6
        for row, column in ((0, 0), (5, 1), (11, 2)):
            perturbed = coordinates.copy(); perturbed[row, column] += eps
            plus = J_of(perturbed)
            perturbed[row, column] -= 2 * eps
            minus = J_of(perturbed)
            fd = (plus - minus) / (2 * eps)
            self.assertAlmostEqual(result["coordinates"][row, column], fd,
                                   delta=1e-5 * max(1.0, abs(fd)))

    def test_NoGradPathUnchanged(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import point_cloud_inference_process
        model = torch.nn.Linear(4, 1).double()
        prediction, _ = point_cloud_inference_process.RunPointCloudForward(
            model, "cpu", "generic",
            torch.rand(5, 1, dtype=torch.float64), torch.rand(5, 3, dtype=torch.float64))
        self.assertFalse(prediction.requires_grad)

    def test_WrtValidation(self):
        model = torch.nn.Linear(4, 1).double()
        with self.assertRaisesRegex(ValueError, "wrt"):
            sensitivity_utils.ComputeSurrogateSensitivities(
                model, numpy.zeros((3, 1)), numpy.zeros((3, 3)),
                lambda p: p.sum(), wrt=("shape",))
        with self.assertRaisesRegex(ValueError, "wrt"):
            sensitivity_utils.ComputeSurrogateSensitivities(
                model, numpy.zeros((3, 1)), numpy.zeros((3, 3)),
                lambda p: p.sum(), wrt=())


@KratosUnittest.skipUnless(have_convection_diffusion,
                           "Missing ConvectionDiffusion/LinearSolvers applications.")
class TestShapeSensitivityFieldOnAnotherPhysics(KratosUnittest.TestCase):
    """The element-local shape field on a DIFFERENT element and physics.

    test_adjoint_cross_validation pins it against StructuralMechanics'
    3D4N solid elements. This guards against the routine being accidentally
    specialized to that case: ConvectionDiffusion is a scalar field, a
    different element type, and a larger mesh. Needs no torch.
    """

    def setUp(self):
        self.model, self.model_part = _SolveThermal(conductivity=2.0, heat_flux=1.0,
                                                    divisions=6)
        self.assembler = differentiable_residual.TangentAssembler(self.model_part)
        self.dof_map = differentiable_residual.DofFieldMap(
            self.assembler, [("TEMPERATURE", "node_historical")])
        # J = sum of nodal temperatures
        self.dJ_du = numpy.ones(self.dof_map.n_equations)

    def test_FieldMatchesThePerCoordinateGlobalPathOnALargerMesh(self):
        field = sensitivity_utils.ComputeShapeSensitivityField(
            self.assembler, self.dof_map, self.dJ_du, fd_step=1e-6)
        self.assertEqual(field.shape, (self.model_part.NumberOfNodes(), 3))
        self.assertGreater(self.model_part.NumberOfElements(), 50)

        # a handful of interior nodes: the global path costs two full
        # re-assemblies per coordinate, so only a sample is affordable - which
        # is the whole reason the element-local pass exists
        rows = {node.Id: row for row, node in enumerate(self.model_part.Nodes)}
        interior = [node for node in self.model_part.Nodes
                    if 1e-9 < node.X0 < 1.0 - 1e-9 and 1e-9 < node.Y0 < 1.0 - 1e-9][:3]
        self.assertEqual(len(interior), 3)

        axes = (("X0", "X"), ("Y0", "Y"), ("Z0", "Z"))
        appliers = {}
        for node in interior:
            for axis in range(2):   # planar case: z is degenerate
                def MakeApplier(node=node, axis=axis):
                    def ApplyCoordinate(value):
                        for attribute in axes[axis]:
                            setattr(node, attribute, value)
                    return ApplyCoordinate
                appliers[(node.Id, axis)] = (MakeApplier(),
                                             getattr(node, axes[axis][0]))
        reference = sensitivity_utils.ComputeParameterSensitivities(
            self.assembler, self.dof_map, self.dJ_du, appliers, fd_step=1e-6)

        for (node_id, axis), value in reference.items():
            with self.subTest(node=node_id, axis=axis):
                self.assertGreater(abs(value), 1e-12)
                self.assertAlmostEqual(field[rows[node_id], axis], value,
                                       delta=1e-6 * abs(value))

    def test_DesignRestrictionCostsLessAndAgrees(self):
        design_ids = [node.Id for node in self.model_part.Nodes
                      if abs(node.X0 - 1.0) < 1e-9]
        self.assertGreater(len(design_ids), 2)

        full = sensitivity_utils.ComputeShapeSensitivityField(
            self.assembler, self.dof_map, self.dJ_du, fd_step=1e-6)
        restricted = sensitivity_utils.ComputeShapeSensitivityField(
            self.assembler, self.dof_map, self.dJ_du, fd_step=1e-6,
            design_node_ids=design_ids)

        rows = {node.Id: row for row, node in enumerate(self.model_part.Nodes)}
        for node_id in design_ids:
            row = rows[node_id]
            for axis in range(3):
                self.assertAlmostEqual(restricted[row, axis], full[row, axis],
                                       delta=1e-9 * abs(full[row, axis]) + 1e-18)
        untouched = [rows[node.Id] for node in self.model_part.Nodes
                     if node.Id not in set(design_ids)]
        self.assertEqual(numpy.abs(restricted[untouched]).max(), 0.0)


if __name__ == '__main__':
    KratosUnittest.main()
