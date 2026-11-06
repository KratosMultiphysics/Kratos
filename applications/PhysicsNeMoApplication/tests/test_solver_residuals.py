"""Tests for solver_residuals: PDE residuals of predicted fields, assembled
by the real solver machinery (ConvectionDiffusion-gated)."""

import sys
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.physics import solver_residuals
have_convection_diffusion = kratos_utils.CheckIfApplicationsAvailable(
    "ConvectionDiffusionApplication", "LinearSolversApplication")

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))


@KratosUnittest.skipUnless(have_convection_diffusion,
                           "Missing required applications: ConvectionDiffusionApplication, LinearSolversApplication.")
class TestResidualEvaluator(KratosUnittest.TestCase):
    def _SolveThermalCase(self, model):
        import thermal_case
        analysis = thermal_case.CreateThermalAnalysis(
            model, conductivity=1.0, heat_flux=1.0, divisions=6)
        analysis.Run()
        return model["ThermalModelPart"]

    def test_ResidualVanishesAtTheSolution(self):
        model = Kratos.Model()
        model_part = self._SolveThermalCase(model)
        evaluator = solver_residuals.BuildResidualEvaluator(model_part)
        self.assertLess(evaluator.ComputeResidualNorm(), 1e-8)

    def test_ResidualGrowsMonotonicallyWithPerturbation(self):
        model = Kratos.Model()
        model_part = self._SolveThermalCase(model)
        evaluator = solver_residuals.BuildResidualEvaluator(model_part)

        solution = {node.Id: node.GetSolutionStepValue(Kratos.TEMPERATURE)
                    for node in model_part.Nodes}

        def perturbed_norm(delta):
            for node in model_part.Nodes:
                if not node.IsFixed(Kratos.TEMPERATURE):
                    node.SetSolutionStepValue(Kratos.TEMPERATURE, solution[node.Id] + delta)
            return evaluator.ComputeResidualNorm()

        base = perturbed_norm(0.0)
        small = perturbed_norm(0.01)
        large = perturbed_norm(0.02)
        self.assertLess(base, 1e-8)
        self.assertGreater(small, 1e-6)
        # linear problem: the residual is exactly linear in the perturbation
        self.assertGreater(large, 1.5 * small)

    def test_NodalResidualsRespectFixedDofs(self):
        model = Kratos.Model()
        model_part = self._SolveThermalCase(model)
        evaluator = solver_residuals.BuildResidualEvaluator(model_part)

        for node in model_part.Nodes:
            if not node.IsFixed(Kratos.TEMPERATURE):
                node.SetSolutionStepValue(
                    Kratos.TEMPERATURE, node.GetSolutionStepValue(Kratos.TEMPERATURE) + 0.5)

        nodal = evaluator.ComputeNodalResiduals()
        self.assertEqual(len(nodal), model_part.NumberOfNodes())
        free_values = []
        for node in model_part.Nodes:
            value = nodal[(node.Id, "TEMPERATURE")]
            if node.IsFixed(Kratos.TEMPERATURE):
                self.assertEqual(value, 0.0)
            else:
                free_values.append(value)
        self.assertTrue(any(value > 1e-8 for value in free_values))


if __name__ == '__main__':
    KratosUnittest.main()
