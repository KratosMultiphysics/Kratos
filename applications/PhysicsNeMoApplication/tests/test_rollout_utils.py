import numpy

import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.training import rollout_utils
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import physicsnemo.metrics.general.mse  # noqa: F401
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestEvaluateRollout(KratosUnittest.TestCase):
    def setUp(self):
        class LinearExtrapolator(torch.nn.Module):
            def forward(self, window):
                width = window.shape[1] // 2
                return 2.0 * window[:, width:] - window[:, :width]

        self.model = torch.jit.script(LinearExtrapolator())

    def _Trajectory(self, fn, steps=10, nodes=5):
        # states (T, N, 1): per-node trajectories u_i(t) = fn(i, t)
        return numpy.array([[[fn(i, t)] for i in range(nodes)] for t in range(steps)])

    def test_ExactOnLinearTrajectory(self):
        states = self._Trajectory(lambda i, t: (1.0 + i) * t)
        predictions, errors = rollout_utils.EvaluateRollout(self.model, states, history_size=2)
        self.assertEqual(predictions.shape, (8, 5, 1))
        self.assertEqual(errors.shape, (8,))
        # Linear extrapolation is exact for linear-in-time fields, even fed
        # its own predictions: zero error at every rollout step.
        self.assertTrue(numpy.allclose(errors, 0.0, atol=1e-12))

    def test_ErrorGrowsOnQuadraticTrajectory(self):
        states = self._Trajectory(lambda i, t: t * t)
        _, errors = rollout_utils.EvaluateRollout(self.model, states, history_size=2)
        # Linear extrapolation of a quadratic drifts, and autoregressive
        # feedback compounds it: strictly growing error curve.
        self.assertGreater(errors[0], 0.0)
        self.assertTrue(numpy.all(numpy.diff(errors) > 0.0))

    @KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
    def test_PerStepMetricsExtras(self):
        states = self._Trajectory(lambda i, t: (1.0 + i) * t)
        result = rollout_utils.EvaluateRollout(
            self.model, states, history_size=2, metric_names=["rmse", "relative_l2"])
        self.assertEqual(len(result), 3)
        predictions, errors, extras = result
        per_step = extras["per_step_metrics"]
        self.assertEqual(per_step["rmse"].shape, errors.shape)
        # the per-step rmse metric must agree with the built-in error curve
        self.assertTrue(numpy.allclose(per_step["rmse"], errors, atol=1e-12))
        self.assertEqual(per_step["relative_l2"].shape, errors.shape)

    def test_MonteCarloRolloutStd(self):
        class DropoutExtrapolator(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.dropout = torch.nn.Dropout(p=0.3)

            def forward(self, window):
                width = window.shape[1] // 2
                return self.dropout(2.0 * window[:, width:] - window[:, :width])

        model = DropoutExtrapolator().double()
        model.eval()
        states = self._Trajectory(lambda i, t: 1.0 + t, steps=6)
        predictions, errors, extras = rollout_utils.EvaluateRollout(
            model, states, history_size=2, mc_samples=8, seed=0)
        self.assertEqual(extras["std"].shape, predictions.shape)
        self.assertGreater(float(numpy.abs(extras["std"]).max()), 0.0)

    def test_InvalidInputsRejected(self):
        states = self._Trajectory(lambda i, t: t)
        with self.assertRaisesRegex(ValueError, "history_size"):
            rollout_utils.EvaluateRollout(self.model, states, history_size=1)
        with self.assertRaisesRegex(ValueError, "shape"):
            rollout_utils.EvaluateRollout(self.model, states[:, :, 0], history_size=2)
        with self.assertRaisesRegex(ValueError, "roll out"):
            rollout_utils.EvaluateRollout(self.model, states[:2], history_size=2)


if __name__ == '__main__':
    KratosUnittest.main()
