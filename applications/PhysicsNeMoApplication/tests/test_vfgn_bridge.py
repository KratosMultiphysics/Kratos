"""Tests for the Virtual Foundry GraphNet (sintering) bridge."""

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges import particle_bridge
from KratosMultiphysics.PhysicsNeMoApplication.training import torch_dataset
from KratosMultiphysics.PhysicsNeMoApplication.bridges import vfgn_bridge
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import torch_scatter  # noqa: F401
    from physicsnemo.models.vfgn import VFGNLearnedSimulator
    have_vfgn = have_torch
except ImportError:
    have_vfgn = False

_MISSING = "Missing required python modules: torch, torch_scatter, physicsnemo."


def _SinteringTrajectory(steps=8, n=6, seed=0):
    """A shrinking particle cloud - sintering's defining behaviour."""
    rng = numpy.random.default_rng(seed)
    initial = rng.random((n, 3)) * 0.1
    centre = initial.mean(axis=0)
    return numpy.stack([centre + (initial - centre) * (1.0 - 0.02 * t)
                        for t in range(steps)])


def _Stats():
    return vfgn_bridge.MakeNormalizationStats(
        velocity_mean=[0.0] * 3, velocity_std=[1.0] * 3,
        acceleration_mean=[0.0] * 3, acceleration_std=[1.0] * 3)


def _Edges(positions, radius=0.05):
    edge_index, _ = particle_bridge.BuildParticleGraphFromPositions(
        positions, Kratos.Parameters(
            '{"type": "radius", "radius": %f, "backend": "numpy"}' % radius))
    return edge_index


class TestNormalizationStats(KratosUnittest.TestCase):
    """Pure numpy/torch: the card round trip needs no physicsnemo."""

    @KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
    def test_CardRoundTrip(self):
        stats = vfgn_bridge.MakeNormalizationStats(
            [0.1, 0.2, 0.3], [1.0, 2.0, 3.0], [0.0] * 3, [0.5] * 3)
        card = vfgn_bridge.StatsToCard(stats)
        restored = vfgn_bridge.StatsFromCard(card)
        numpy.testing.assert_allclose(
            numpy.asarray(restored["velocity"].mean), [0.1, 0.2, 0.3], atol=1e-12)
        numpy.testing.assert_allclose(
            numpy.asarray(restored["acceleration"].std), [0.5] * 3, atol=1e-12)

    @KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
    def test_ZeroStdRejected(self):
        """Upstream divides by std without an epsilon for velocity and
        acceleration, so a zero would become a silent NaN."""
        with self.assertRaisesRegex(ValueError, "standard deviation"):
            vfgn_bridge.MakeNormalizationStats(
                [0.0] * 3, [1.0, 0.0, 1.0], [0.0] * 3, [1.0] * 3)

    @KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
    def test_StatsFromTheShippedTrajectoryDataset(self):
        trajectory = _SinteringTrajectory(steps=10)
        dataset = torch_dataset.CreateParticleTrajectoryDataset(
            trajectory, history_size=2, delta_time=0.1)
        stats = vfgn_bridge.StatsFromTrajectoryDataset(dataset)
        self.assertEqual(tuple(numpy.asarray(stats["velocity"].mean).shape), (3,))
        self.assertEqual(tuple(numpy.asarray(stats["acceleration"].std).shape), (3,))

    @KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
    def test_CardWithoutStatsRejected(self):
        with self.assertRaisesRegex(ValueError, "normalization statistics"):
            vfgn_bridge.StatsFromCard({})


@KratosUnittest.skipUnless(have_vfgn, _MISSING)
class TestVfgnSimulator(KratosUnittest.TestCase):

    def setUp(self):
        self.trajectory = _SinteringTrajectory()
        self.settings = Kratos.Parameters("""{
            "predict_length"     : 2,
            "num_seq"            : 5,
            "num_particle_types" : 3,
            "connectivity_param" : 0.05,
            "boundaries"         : [[0.0, 0.2], [0.0, 0.2], [0.0, 0.2]]
        }""")

    def test_RolloutProducesFinitePositions(self):
        torch.manual_seed(0)
        model = vfgn_bridge.CreateVfgnSimulator(self.settings, _Stats())
        sequence = self.trajectory[:5]                       # (T, N, 3)
        sequence = numpy.transpose(sequence, (1, 0, 2))      # -> (N, T, 3)
        edges = _Edges(sequence[:, -1, :])

        predicted = vfgn_bridge.RunVfgnRollout(model, sequence, edges, predict_length=2)
        self.assertEqual(predicted.shape, (sequence.shape[0], 2, 3))
        self.assertTrue(numpy.isfinite(predicted).all())

    def test_TrainingStepReducesTheLoss(self):
        torch.manual_seed(0)
        model = vfgn_bridge.CreateVfgnSimulator(self.settings, _Stats())
        sequence = numpy.transpose(self.trajectory[:5], (1, 0, 2))
        next_positions = numpy.transpose(self.trajectory[5:7], (1, 0, 2))
        edges = _Edges(sequence[:, -1, :])

        # the optimizer must come AFTER a forward: the model creates
        # parameters lazily and would otherwise optimize about half of them
        first = vfgn_bridge.ComputeVfgnLoss(
            model, sequence, next_positions, edges, predict_length=2)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

        history = [float(first)]
        for _ in range(12):
            optimizer.zero_grad()
            loss = vfgn_bridge.ComputeVfgnLoss(
                model, sequence, next_positions, edges, predict_length=2)
            loss.backward()
            optimizer.step()
            history.append(float(loss))
        # the input noise is redrawn every call, so compare against the best
        self.assertLess(min(history[1:]), history[0])

    def test_LazyParametersAppearOnFirstForward(self):
        torch.manual_seed(0)
        model = vfgn_bridge.CreateVfgnSimulator(self.settings, _Stats())
        before = sum(p.numel() for p in model.parameters())
        sequence = numpy.transpose(self.trajectory[:5], (1, 0, 2))
        vfgn_bridge.RunVfgnRollout(model, sequence, _Edges(sequence[:, -1, :]),
                                   predict_length=2)
        after = sum(p.numel() for p in model.parameters())
        self.assertGreater(after, before)

    def test_Validation(self):
        with self.assertRaisesRegex(ValueError, "predict_length"):
            vfgn_bridge.CreateVfgnSimulator(
                Kratos.Parameters('{"predict_length": 0}'), _Stats())
        with self.assertRaisesRegex(ValueError, "num_seq"):
            vfgn_bridge.CreateVfgnSimulator(
                Kratos.Parameters('{"num_seq": 2}'), _Stats())
        with self.assertRaisesRegex(ValueError, "normalization statistics"):
            vfgn_bridge.CreateVfgnSimulator(self.settings, None)

        model = vfgn_bridge.CreateVfgnSimulator(self.settings, _Stats())
        sequence = numpy.transpose(self.trajectory[:5], (1, 0, 2))
        with self.assertRaisesRegex(ValueError, r"\(N, T, 3\)"):
            vfgn_bridge.RunVfgnRollout(model, sequence[:, :, :2],
                                       _Edges(sequence[:, -1, :]), predict_length=2)
        with self.assertRaisesRegex(ValueError, "predict_length"):
            vfgn_bridge.ComputeVfgnLoss(
                model, sequence, numpy.transpose(self.trajectory[5:6], (1, 0, 2)),
                _Edges(sequence[:, -1, :]), predict_length=2)


@KratosUnittest.skipUnless(have_vfgn, _MISSING)
class TestUpstreamForwardIsBroken(KratosUnittest.TestCase):
    """Pins the physicsnemo 2.2 bug this bridge routes around.

    VFGNLearnedSimulator.forward()'s shape guard demands a 2-D
    next_positions while its body needs (N, predict_length, 3), so no input
    works. When this test starts failing, upstream has fixed it and
    ComputeVfgnLoss can call forward() directly.
    """

    def test_ForwardRejectsBothShapes(self):
        torch.manual_seed(0)
        model = vfgn_bridge.CreateVfgnSimulator(Kratos.Parameters("""{
            "predict_length": 2, "num_seq": 5, "connectivity_param": 0.05
        }"""), _Stats())
        trajectory = _SinteringTrajectory()
        sequence = torch.as_tensor(numpy.transpose(trajectory[:5], (1, 0, 2)))
        n = sequence.shape[0]
        edges = torch.as_tensor(_Edges(sequence[:, -1, :].numpy()))
        noise = model.get_random_walk_noise_for_position_sequence(sequence, 6.7e-4)
        arguments = (noise, sequence, torch.tensor([n]),
                     torch.tensor([int(edges.shape[1])]), edges[0], edges[1], 2,
                     None, torch.zeros(n, dtype=torch.int64))

        # the guard-legal 2-D shape fails in the body's arithmetic ...
        with self.assertRaises(RuntimeError):
            model.forward(torch.rand(n, 3, dtype=torch.float64), *arguments)
        # ... and the math-legal 3-D shape is rejected by the guard
        with self.assertRaisesRegex(ValueError, "2D"):
            model.forward(torch.rand(n, 2, 3, dtype=torch.float64), *arguments)


if __name__ == '__main__':
    KratosUnittest.main()
