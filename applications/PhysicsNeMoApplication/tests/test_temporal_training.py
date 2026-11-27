"""Tests for the temporal training schemes (window datasets, BPTT rollout
training) on synthetic and real Kratos transient trajectories."""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.training import temporal_training
sys.path.insert(0, str(Path(__file__).parent / "kratos_solver_cases"))

have_structural = kratos_utils.CheckIfApplicationsAvailable(
    "StructuralMechanicsApplication", "LinearSolversApplication")

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False


def _DecayTrajectory(steps=12, nodes=5, width=2, rate=0.8):
    """u_{t+1} = rate * u_t: a linear system a small model can learn exactly."""
    rng = numpy.random.default_rng(0)
    state = rng.standard_normal((nodes, width))
    states = [state]
    for _ in range(steps - 1):
        states.append(rate * states[-1])
    return numpy.stack(states)


def _MakeModel(input_width, output_width, seed=0, hidden=16):
    torch.manual_seed(seed)
    return torch.nn.Sequential(
        torch.nn.Linear(input_width, hidden), torch.nn.Tanh(),
        torch.nn.Linear(hidden, output_width))


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestTrajectoryWindowDataset(KratosUnittest.TestCase):

    def setUp(self):
        self.states = _DecayTrajectory(steps=6, nodes=3, width=2)

    def test_SingleStepWindowsMatchHandBuilt(self):
        dataset = temporal_training.CreateTrajectoryWindowDataset(
            self.states, Kratos.Parameters('{"scheme": "single_step", "history_size": 2}'))
        self.assertEqual(len(dataset), 4)  # T - K
        inputs, targets = dataset[0]
        self.assertEqual(tuple(inputs.shape), (3, 4))   # (N, K*W)
        self.assertEqual(tuple(targets.shape), (3, 2))
        # oldest first: [state_0 | state_1] -> state_2
        expected = numpy.concatenate([self.states[0], self.states[1]], axis=1)
        numpy.testing.assert_allclose(inputs.numpy(), expected, rtol=1e-6)
        numpy.testing.assert_allclose(targets.numpy(), self.states[2], rtol=1e-6)

    def test_OneShotAndTimeConditional(self):
        one_shot = temporal_training.CreateTrajectoryWindowDataset(
            self.states, Kratos.Parameters('{"scheme": "one_shot", "history_size": 2}'))
        self.assertEqual(len(one_shot), 1)
        inputs, targets = one_shot[0]
        numpy.testing.assert_allclose(targets.numpy(), self.states[-1], rtol=1e-6)

        conditional = temporal_training.CreateTrajectoryWindowDataset(
            self.states, Kratos.Parameters('{"scheme": "time_conditional", "history_size": 2}'))
        self.assertEqual(len(conditional), 4)
        inputs, targets = conditional[-1]
        self.assertEqual(tuple(inputs.shape), (3, 5))  # K*W + the time channel
        self.assertAlmostEqual(float(inputs[0, -1]), 1.0, places=6)  # last step -> t = 1
        numpy.testing.assert_allclose(targets.numpy(), self.states[-1], rtol=1e-6)

    def test_MultipleTrajectoriesConcatenate(self):
        dataset = temporal_training.CreateTrajectoryWindowDataset(
            [self.states, self.states],
            Kratos.Parameters('{"scheme": "single_step", "history_size": 2}'))
        self.assertEqual(len(dataset), 8)

    def test_Validation(self):
        with self.assertRaisesRegex(ValueError, "Unknown scheme"):
            temporal_training.CreateTrajectoryWindowDataset(
                self.states, Kratos.Parameters('{"scheme": "teacher_forced"}'))
        with self.assertRaisesRegex(ValueError, "history_size"):
            temporal_training.CreateTrajectoryWindowDataset(
                self.states, Kratos.Parameters('{"history_size": 0}'))
        with self.assertRaisesRegex(ValueError, "history_size"):
            temporal_training.CreateTrajectoryWindowDataset(
                self.states[:2], Kratos.Parameters('{"history_size": 5}'))
        with self.assertRaisesRegex(ValueError, r"\(T, N, W\)"):
            temporal_training.CreateTrajectoryWindowDataset(
                [numpy.zeros((3, 4))], Kratos.Parameters("{}"))


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestAutoregressiveTraining(KratosUnittest.TestCase):

    def setUp(self):
        self.states = _DecayTrajectory(steps=12, nodes=5, width=2)

    def test_CheckpointedGradientsMatchPlain(self):
        """The BPTT correctness pin: checkpointing must not change gradients."""
        seed_states = [torch.tensor(self.states[i], dtype=torch.float32) for i in range(2)]
        targets = torch.tensor(self.states[2:6], dtype=torch.float32)

        gradients = []
        for checkpoint in (False, True):
            model = _MakeModel(4, 2, seed=0)
            model.train()
            predictions = temporal_training.RolloutPredictions(
                model, seed_states, steps=4, checkpoint=checkpoint)
            torch.nn.functional.mse_loss(predictions, targets).backward()
            gradients.append([p.grad.clone() for p in model.parameters()])

        for plain, checkpointed in zip(*gradients):
            self.assertTrue(torch.allclose(plain, checkpointed, atol=1e-6))

    def test_RolloutShapeAndGraph(self):
        model = _MakeModel(4, 2)
        model.train()
        seed_states = [torch.tensor(self.states[i], dtype=torch.float32) for i in range(2)]
        predictions = temporal_training.RolloutPredictions(model, seed_states, steps=3)
        self.assertEqual(tuple(predictions.shape), (3, 5, 2))
        self.assertTrue(predictions.requires_grad)

    def test_TrainingReducesRolloutError(self):
        model = _MakeModel(4, 2, seed=1)
        history = temporal_training.TrainAutoregressive(model, self.states, Kratos.Parameters("""{
            "epochs"        : 60,
            "rollout_steps" : 3,
            "history_size"  : 2,
            "learning_rate" : 5e-3,
            "device"        : "cpu",
            "seed"          : 0
        }"""))
        self.assertEqual(len(history), 60)
        self.assertLess(history[-1], history[0])

        from KratosMultiphysics.PhysicsNeMoApplication.training import rollout_utils
        _, errors = rollout_utils.EvaluateRollout(model, self.states, history_size=2)
        untrained = _MakeModel(4, 2, seed=2)
        _, baseline = rollout_utils.EvaluateRollout(untrained, self.states, history_size=2)
        self.assertLess(float(numpy.mean(errors)), float(numpy.mean(baseline)))

    def test_CheckpointedTrainingRuns(self):
        model = _MakeModel(4, 2, seed=1)
        history = temporal_training.TrainAutoregressive(model, self.states, Kratos.Parameters("""{
            "epochs"                 : 5,
            "rollout_steps"          : 3,
            "gradient_checkpointing" : true,
            "device"                 : "cpu",
            "seed"                   : 0
        }"""))
        self.assertEqual(len(history), 5)
        self.assertTrue(all(numpy.isfinite(history)))

    def test_Validation(self):
        model = _MakeModel(4, 2)
        with self.assertRaisesRegex(ValueError, "rollout_steps"):
            temporal_training.TrainAutoregressive(
                model, self.states, Kratos.Parameters('{"rollout_steps": 0}'))
        with self.assertRaisesRegex(ValueError, "Unknown optimizer"):
            temporal_training.TrainAutoregressive(
                model, self.states, Kratos.Parameters('{"optimizer": "lbfgs", "epochs": 1}'))
        with self.assertRaisesRegex(ValueError, "rollout_steps"):
            temporal_training.TrainAutoregressive(
                model, self.states[:3], Kratos.Parameters('{"rollout_steps": 8, "epochs": 1}'))


@KratosUnittest.skipUnless(have_torch and have_structural,
                           "Missing torch or StructuralMechanicsApplication.")
class TestTransientStructuralSurrogate(KratosUnittest.TestCase):

    def test_SingleStepSurrogateOnRealTransientSolve(self):
        import structural_case
        import transient_harness
        from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
        from KratosMultiphysics.PhysicsNeMoApplication.training import rollout_utils
        model = Kratos.Model()
        analysis = structural_case.CreateTransientStructuralAnalysis(
            model, divisions=4, time_step=0.005, end_time=0.06)
        states = transient_harness.RunTransientAnalysis(
            analysis, collect=structural_case.CollectDisplacements)
        self.assertEqual(states.ndim, 3)
        self.assertGreaterEqual(states.shape[0], 10)
        # a genuine oscillation, not a monotone ramp
        tip = states[:, -1, 1]
        self.assertGreater(float(numpy.std(tip)), 0.0)

        # normalize: raw displacements are ~1e-3, too small to train on directly
        scale = float(numpy.abs(states).max())
        scaled = states / scale

        dataset = temporal_training.CreateTrajectoryWindowDataset(
            scaled, Kratos.Parameters('{"scheme": "single_step", "history_size": 2}'))
        self.assertEqual(len(dataset), states.shape[0] - 2)

        width = states.shape[2]
        surrogate = _MakeModel(2 * width, width, seed=3, hidden=32)
        history = training_utils.TrainModel(surrogate, dataset, Kratos.Parameters("""{
            "epochs"        : 250,
            "batch_size"    : 4,
            "learning_rate" : 5e-3,
            "device"        : "cpu",
            "seed"          : 0
        }"""))
        self.assertLess(history[-1], history[0])

        _, errors = rollout_utils.EvaluateRollout(surrogate, scaled, history_size=2)
        untrained = _MakeModel(2 * width, width, seed=4, hidden=32)
        _, baseline = rollout_utils.EvaluateRollout(untrained, scaled, history_size=2)
        self.assertLess(float(numpy.mean(errors)), float(numpy.mean(baseline)))


if __name__ == '__main__':
    KratosUnittest.main()
