"""Tests for differentiable_residual: the assembled FEM residual as a
torch.autograd.Function (ConvectionDiffusion/StructuralMechanics-gated)."""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.physics import differentiable_residual
from KratosMultiphysics.PhysicsNeMoApplication.physics import solver_residuals
from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

have_convection_diffusion = kratos_utils.CheckIfApplicationsAvailable(
    "ConvectionDiffusionApplication", "LinearSolversApplication")
have_structural = kratos_utils.CheckIfApplicationsAvailable(
    "StructuralMechanicsApplication", "LinearSolversApplication")

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))

_THERMAL_FIELDS = [("TEMPERATURE", "node_historical")]


def _SolveThermalCase(model, divisions=3):
    import thermal_case
    analysis = thermal_case.CreateThermalAnalysis(
        model, conductivity=1.0, heat_flux=1.0, divisions=divisions)
    analysis.Run()
    return model["ThermalModelPart"]


@KratosUnittest.skipUnless(have_torch and have_convection_diffusion,
                           "Missing torch or ConvectionDiffusion/LinearSolvers applications.")
class TestDifferentiableResidual(KratosUnittest.TestCase):
    def _Setup(self, divisions=3):
        model = Kratos.Model()
        model_part = _SolveThermalCase(model, divisions)
        assembler = differentiable_residual.TangentAssembler(model_part)
        dof_map = differentiable_residual.DofFieldMap(assembler, _THERMAL_FIELDS)
        return model, model_part, assembler, dof_map

    def test_DofFieldMapRoundTrip(self):
        _, model_part, assembler, dof_map = self._Setup()
        u = dof_map.ReadDofVector()
        self.assertEqual(u.shape, (dof_map.n_equations,))
        # fields -> dofs -> fields is the identity on the DOF entries
        fields = dof_map.DofVectorToFields(u)
        self.assertEqual(fields.shape, (model_part.NumberOfNodes(), 1))
        numpy.testing.assert_array_equal(dof_map.FieldsToDofVector(fields), u)
        # write/read round trip
        perturbed = u + 0.123
        dof_map.WriteDofVector(perturbed)
        numpy.testing.assert_allclose(dof_map.ReadDofVector(), perturbed, atol=1e-14)
        # the nodal database saw the write
        temperatures = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in model_part.Nodes])
        numpy.testing.assert_allclose(
            numpy.sort(temperatures), numpy.sort(perturbed), atol=1e-14)

    def test_ResidualMatchesPlainEvaluator(self):
        _, _, assembler, dof_map = self._Setup()
        u0 = torch.tensor(dof_map.ReadDofVector(), dtype=torch.float64)
        b = differentiable_residual.KratosResidualFunction.Apply(u0, assembler, dof_map)
        self.assertLess(float(b.abs().max()), 1e-8)  # at the solution: b = 0

    def test_SignConventionPinned(self):
        # linear problem: b(u0 + delta) - b(u0) == -(masked K) delta, exactly
        _, _, assembler, dof_map = self._Setup()
        u0 = dof_map.ReadDofVector()

        dof_map.WriteDofVector(u0)
        b0 = numpy.array(assembler.ComputeResidualVector(), copy=True)
        K = assembler.ComputeTangentMatrix(apply_dirichlet=False)

        rng = numpy.random.default_rng(0)
        delta = rng.standard_normal(u0.shape)
        dof_map.WriteDofVector(u0 + delta)
        b1 = numpy.array(assembler.ComputeResidualVector(), copy=True)

        predicted = -(K @ delta)
        predicted[dof_map.fixed_mask] = 0.0  # BuildRHS zeroes fixed rows
        numpy.testing.assert_allclose(b1 - b0, predicted, atol=1e-9)

        # and backward of b . g must be -(K^T (masked g))
        u_t = torch.tensor(u0, dtype=torch.float64, requires_grad=True)
        b = differentiable_residual.KratosResidualFunction.Apply(u_t, assembler, dof_map)
        g = torch.tensor(rng.standard_normal(u0.shape), dtype=torch.float64)
        (b * g).sum().backward()
        g_masked = g.numpy().copy()
        g_masked[dof_map.fixed_mask] = 0.0
        numpy.testing.assert_allclose(u_t.grad.numpy(), -(K.T @ g_masked), atol=1e-9)

    def test_Gradcheck(self):
        _, _, assembler, dof_map = self._Setup(divisions=2)
        u0 = torch.tensor(dof_map.ReadDofVector(), dtype=torch.float64, requires_grad=True)
        self.assertTrue(torch.autograd.gradcheck(
            lambda u: differentiable_residual.KratosResidualFunction.Apply(
                u, assembler, dof_map),
            (u0,), eps=1e-6, atol=1e-7, rtol=1e-5))

    def test_FdDirectionalDerivative(self):
        _, _, assembler, dof_map = self._Setup()
        rng = numpy.random.default_rng(1)
        u0 = dof_map.ReadDofVector() + 0.05 * rng.standard_normal(dof_map.n_equations)
        direction = rng.standard_normal(dof_map.n_equations)

        def loss_of(u_numpy):
            u = torch.tensor(u_numpy, dtype=torch.float64, requires_grad=True)
            b = differentiable_residual.KratosResidualFunction.Apply(u, assembler, dof_map)
            return b.square().mean(), u

        loss, u = loss_of(u0)
        loss.backward()
        autograd_directional = float(u.grad.numpy() @ direction)

        eps = 1e-6
        loss_plus, _ = loss_of(u0 + eps * direction)
        loss_minus, _ = loss_of(u0 - eps * direction)
        fd_directional = (float(loss_plus) - float(loss_minus)) / (2.0 * eps)
        self.assertAlmostEqual(autograd_directional, fd_directional,
                               delta=1e-6 * max(1.0, abs(fd_directional)))

    def test_ExactResidualLossTrainsTinyModel(self):
        model, model_part, _, _ = self._Setup(divisions=3)

        # inputs: node coordinates; targets: a deliberately WRONG field
        # (scaled solution) - the exact-residual term must pull the model
        # toward the true physics despite the biased data
        coordinates = torch.tensor(
            [[node.X, node.Y] for node in model_part.Nodes], dtype=torch.float64)
        solution = torch.tensor(
            [[node.GetSolutionStepValue(Kratos.TEMPERATURE)] for node in model_part.Nodes],
            dtype=torch.float64)
        dataset = torch.utils.data.TensorDataset(coordinates, 0.5 * solution)

        torch.manual_seed(0)
        surrogate = torch.nn.Sequential(
            torch.nn.Linear(2, 32), torch.nn.Tanh(), torch.nn.Linear(32, 1)).double()

        term = differentiable_residual.MakeExactResidualLossTerm(
            Kratos.Parameters("""{
                "fields" : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "weight" : 5.0
            }"""),
            model_part, lambda: coordinates)

        evaluator = solver_residuals.BuildResidualEvaluator(model_part)

        def residual_of(model_out):
            for node, value in zip(model_part.Nodes, model_out.reshape(-1).tolist()):
                node.SetSolutionStepValue(Kratos.TEMPERATURE, value)
            return evaluator.ComputeResidualNorm()

        with torch.no_grad():
            initial = residual_of(surrogate(coordinates))
        training_utils.TrainModel(surrogate, dataset, Kratos.Parameters("""{
            "epochs"        : 60,
            "batch_size"    : 64,
            "learning_rate" : 5e-3,
            "device"        : "cpu",
            "seed"          : 0
        }"""), extra_loss_terms=[term])
        with torch.no_grad():
            final = residual_of(surrogate(coordinates))
        self.assertLess(final, initial)


@KratosUnittest.skipUnless(have_torch and have_structural,
                           "Missing torch or StructuralMechanics/LinearSolvers applications.")
class TestDifferentiableResidualStructural(KratosUnittest.TestCase):
    def test_GradcheckMultiComponent(self):
        # pins the DISPLACEMENT_X/Y component mapping of DofFieldMap
        import structural_case
        import KratosMultiphysics.StructuralMechanicsApplication as SMA

        model = Kratos.Model()
        analysis = structural_case.CreateStructuralAnalysis(model, divisions=2)
        analysis.Run()
        model_part = model["StructuralModelPart"]

        assembler = differentiable_residual.TangentAssembler(model_part)
        dof_map = differentiable_residual.DofFieldMap(
            assembler, [("DISPLACEMENT", "node_historical")])
        self.assertEqual(dof_map.total_width, 3)

        u0 = torch.tensor(dof_map.ReadDofVector(), dtype=torch.float64, requires_grad=True)
        b0 = differentiable_residual.KratosResidualFunction.Apply(u0, assembler, dof_map)
        # at the solution the residual is round-off relative to the load scale
        # (the residual at zero displacement = the external forces)
        b_zero = differentiable_residual.KratosResidualFunction.Apply(
            torch.zeros_like(u0), assembler, dof_map)
        self.assertLess(float(b0.detach().abs().max()),
                        1e-6 * float(b_zero.detach().abs().max()))

        # Gradcheck on a UNIT-SCALE material. With the real 210 GPa modulus
        # the tangent entries are ~1e11, and gradcheck's eps=1e-6 finite
        # differences amplify the parallel assembly's summation-order noise
        # by 1/eps - which made this check flaky (~1 run in 3). The mapping
        # under test is scale-free, so rescale instead of loosening
        # tolerances.
        properties = next(iter(model_part.Elements)).Properties
        properties.SetValue(Kratos.YOUNG_MODULUS, 1.0)
        for node in model_part.Nodes:
            node.SetSolutionStepValue(SMA.POINT_LOAD, [0.0, 0.0, 0.0])
        scaled_assembler = differentiable_residual.TangentAssembler(model_part)
        scaled_map = differentiable_residual.DofFieldMap(
            scaled_assembler, [("DISPLACEMENT", "node_historical")])
        u_scaled = torch.tensor(
            numpy.linspace(0.0, 0.1, scaled_map.n_equations),
            dtype=torch.float64, requires_grad=True)
        self.assertTrue(torch.autograd.gradcheck(
            lambda u: differentiable_residual.KratosResidualFunction.Apply(
                u, scaled_assembler, scaled_map),
            (u_scaled,), eps=1e-6, atol=1e-8, rtol=1e-6))


@KratosUnittest.skipUnless(have_torch and have_convection_diffusion,
                           "Missing torch or ConvectionDiffusion/LinearSolvers applications.")
class TestTransientResidualElementIntegrated(KratosUnittest.TestCase):
    """ConvectionDiffusion's transient solver integrates in the ELEMENT (its
    scheme is the static one), so the shipped assembler handles it unchanged
    - the time dependence enters through ProcessInfo and the buffer."""

    def _SolveSomeSteps(self, divisions=3):
        import thermal_case
        import transient_harness

        model = Kratos.Model()
        analysis = thermal_case.CreateTransientThermalAnalysis(
            model, divisions=divisions, time_step=0.05, end_time=0.15)
        states = transient_harness.RunTransientAnalysis(analysis, collect=lambda mp: [
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in mp.Nodes])
        model_part = model["ThermalModelPart"]
        assembler = differentiable_residual.TangentAssembler(model_part)
        dof_map = differentiable_residual.DofFieldMap(assembler, _THERMAL_FIELDS)
        return states, model_part, assembler, dof_map

    def test_TrajectoryIsTransient(self):
        states, _, _, _ = self._SolveSomeSteps()
        self.assertEqual(states.shape[0], 3)
        # the field grows step by step towards the steady state
        self.assertLess(states[0].max(), states[-1].max())

    def test_SignConventionAtMidTrajectoryStep(self):
        _, _, assembler, dof_map = self._SolveSomeSteps()
        assembler.InitializeSolutionStep()
        u0 = dof_map.ReadDofVector()

        dof_map.WriteDofVector(u0)
        b0 = numpy.array(assembler.ComputeResidualVector(), copy=True)
        K = assembler.ComputeTangentMatrix(apply_dirichlet=False)

        rng = numpy.random.default_rng(0)
        delta = rng.standard_normal(u0.shape)
        dof_map.WriteDofVector(u0 + delta)
        b1 = numpy.array(assembler.ComputeResidualVector(), copy=True)

        predicted = -(K @ delta)
        predicted[dof_map.fixed_mask] = 0.0
        numpy.testing.assert_allclose(b1 - b0, predicted, atol=1e-9)

    def test_ResidualIsSmallAtTheSolvedState(self):
        _, _, assembler, dof_map = self._SolveSomeSteps()
        u0 = torch.tensor(dof_map.ReadDofVector(), dtype=torch.float64)
        b_solution = differentiable_residual.KratosResidualFunction.Apply(
            u0, assembler, dof_map)
        b_zero = differentiable_residual.KratosResidualFunction.Apply(
            torch.zeros_like(u0), assembler, dof_map)
        self.assertLess(float(b_solution.abs().max()),
                        1e-8 * float(b_zero.abs().max()))

    def test_Gradcheck(self):
        _, _, assembler, dof_map = self._SolveSomeSteps(divisions=2)
        u0 = torch.tensor(dof_map.ReadDofVector(), dtype=torch.float64, requires_grad=True)
        self.assertTrue(torch.autograd.gradcheck(
            lambda u: differentiable_residual.KratosResidualFunction.Apply(
                u, assembler, dof_map),
            (u0,), eps=1e-6, atol=1e-7, rtol=1e-5))


@KratosUnittest.skipUnless(have_torch and have_structural,
                           "Missing torch or StructuralMechanics/LinearSolvers applications.")
class TestTransientResidualDynamicScheme(KratosUnittest.TestCase):
    """Structural dynamics integrates in the SCHEME (Bossak), so the
    assembler needs scheme= plus the per-step InitializeSolutionStep."""

    def _SolveSomeSteps(self, divisions=3):
        import structural_case
        import transient_harness

        model = Kratos.Model()
        analysis = structural_case.CreateTransientStructuralAnalysis(
            model, divisions=divisions, time_step=0.005, end_time=0.02)
        transient_harness.RunTransientAnalysis(analysis)
        model_part = model["StructuralModelPart"]

        scheme = Kratos.ResidualBasedBossakDisplacementScheme(-0.3)
        assembler = differentiable_residual.TangentAssembler(model_part, scheme=scheme)
        assembler.InitializeSolutionStep()
        dof_map = differentiable_residual.DofFieldMap(
            assembler, [("DISPLACEMENT", "node_historical")])
        return model_part, assembler, dof_map

    def test_SignConventionWithEffectiveTangent(self):
        _, assembler, dof_map = self._SolveSomeSteps()
        u0 = dof_map.ReadDofVector()

        dof_map.WriteDofVector(u0)
        b0 = numpy.array(assembler.ComputeResidualVector(), copy=True)
        K_effective = assembler.ComputeTangentMatrix(apply_dirichlet=False)

        rng = numpy.random.default_rng(0)
        delta = 1e-6 * rng.standard_normal(u0.shape)
        dof_map.WriteDofVector(u0 + delta)
        b1 = numpy.array(assembler.ComputeResidualVector(), copy=True)

        predicted = -(K_effective @ delta)
        predicted[dof_map.fixed_mask] = 0.0
        scale = max(float(numpy.abs(predicted).max()), 1e-30)
        numpy.testing.assert_allclose(b1 - b0, predicted, atol=1e-6 * scale)

    def test_MassTermMakesTheTangentDifferFromStatics(self):
        model_part, assembler, _ = self._SolveSomeSteps()
        K_effective = assembler.ComputeTangentMatrix(apply_dirichlet=False)
        static_assembler = differentiable_residual.TangentAssembler(model_part)
        K_static = static_assembler.ComputeTangentMatrix(apply_dirichlet=False)
        # K_eff = K + M (1-alpha) c0 + D c1 with c0 = 1/(beta dt^2): huge here
        self.assertGreater(abs(K_effective - K_static).max(), 0.0)

    def test_Gradcheck(self):
        _, assembler, dof_map = self._SolveSomeSteps(divisions=2)
        u0 = torch.tensor(dof_map.ReadDofVector(), dtype=torch.float64, requires_grad=True)
        self.assertTrue(torch.autograd.gradcheck(
            lambda u: differentiable_residual.KratosResidualFunction.Apply(
                u, assembler, dof_map),
            (u0,), eps=1e-6, atol=1e-4, rtol=1e-3))


if __name__ == '__main__':
    KratosUnittest.main()
