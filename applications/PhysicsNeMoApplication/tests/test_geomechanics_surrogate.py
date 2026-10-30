"""Consolidation: a coupled geomechanics transient and a surrogate over it.

GeoMechanicsApplication is compiled in the reference build and nothing here
had driven it. Consolidation is the natural case: pore water pressure decays
as water drains, which is a transient field a temporal surrogate predicts.
"""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))

have_geomechanics = kratos_utils.CheckIfApplicationsAvailable(
    "GeoMechanicsApplication", "StructuralMechanicsApplication",
    "LinearSolversApplication")

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False


@KratosUnittest.skipUnless(have_geomechanics, "Missing GeoMechanics applications.")
class TestConsolidationCase(KratosUnittest.TestCase):
    def setUp(self):
        import geomechanics_case

        if not geomechanics_case.IsAvailable():
            self.skipTest("GeoMechanics consolidation data is not present.")

    def test_TheColumnSolvesAndCarriesPorePressure(self):
        import geomechanics_case

        model = Kratos.Model()
        analysis = geomechanics_case.CreateConsolidationAnalysis(model, stage=1)
        analysis.Run()
        model_part = geomechanics_case.GetModelPart(model)
        pressures = geomechanics_case.WaterPressures(model_part)

        self.assertEqual(model_part.NumberOfNodes(), 201)
        self.assertTrue(numpy.isfinite(pressures).all())
        self.assertGreater(numpy.abs(pressures).max(), 0.0)

    def test_TheStagesAreTheTimeAxis(self):
        """This benchmark's transient is a SEQUENCE OF STAGES, each
        continuing the last on the same model - so the series comes from
        running them in order, not from stepping one analysis."""
        import geomechanics_case

        model = Kratos.Model()
        states, model_part = geomechanics_case.SolveStages(model, stages=(1, 2, 3))
        self.assertEqual(states.shape, (3, model_part.NumberOfNodes()))
        self.assertTrue(numpy.isfinite(states).all())
        # consolidation is a decay: the field must actually change
        self.assertGreater(numpy.abs(states[-1] - states[0]).max(), 1e-9)

    def test_TheFixtureIsReusedNotCopied(self):
        import geomechanics_case

        self.assertIn("GeoMechanicsApplication", str(geomechanics_case._GEO_TESTS))


@KratosUnittest.skipUnless(have_geomechanics and have_torch,
                           "Missing GeoMechanics applications or torch.")
class TestConsolidationSurrogate(KratosUnittest.TestCase):
    def setUp(self):
        import geomechanics_case

        if not geomechanics_case.IsAvailable():
            self.skipTest("GeoMechanics consolidation data is not present.")
        self.states, self.model_part = geomechanics_case.SolveStages(
            Kratos.Model(), stages=(1, 2, 3, 4))

    def test_TheSeriesFeedsATemporalDataset(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import temporal_training

        # (T, N, W): one trajectory of a scalar field
        states = torch.tensor(self.states[..., None], dtype=torch.float32)
        dataset = temporal_training.CreateTrajectoryWindowDataset(
            [states], Kratos.Parameters("{}"))
        self.assertGreater(len(dataset), 0)
        inputs, targets = dataset[0]
        self.assertTrue(torch.isfinite(inputs).all())
        self.assertTrue(torch.isfinite(targets).all())

    def test_AOneStepSurrogateBeatsPersistence(self):
        """The honest bar for a transient surrogate: predicting the next
        state must beat repeating the current one."""
        from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils

        states = torch.tensor(self.states, dtype=torch.float32)
        scale = states.abs().max()
        normalized = states / scale
        inputs, targets = normalized[:-1], normalized[1:]

        torch.manual_seed(0)
        model = torch.nn.Sequential(
            torch.nn.Linear(inputs.shape[1], 32), torch.nn.Tanh(),
            torch.nn.Linear(32, inputs.shape[1]))
        history = training_utils.TrainModel(
            model, torch.utils.data.TensorDataset(inputs, targets),
            Kratos.Parameters("""{
                "epochs"        : 300,
                "batch_size"    : 4,
                "learning_rate" : 1e-2,
                "device"        : "cpu",
                "shuffle"       : false,
                "seed"          : 0
            }"""))
        self.assertLess(history[-1], history[0])

        with torch.no_grad():
            predicted = model(inputs)
        surrogate_error = float((predicted - targets).square().mean())
        persistence_error = float((inputs - targets).square().mean())
        Kratos.Logger.PrintInfo(
            "TestConsolidationSurrogate",
            f"one-step MSE {surrogate_error:.3e} against persistence {persistence_error:.3e}")
        self.assertLess(surrogate_error, persistence_error)


if __name__ == '__main__':
    KratosUnittest.main()
