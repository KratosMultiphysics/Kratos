from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.training import rom_temporal
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import physicsnemo.models.mesh_reduced.temporal_model  # noqa: F401
    from physicsnemo.distributed.manager import DistributedManager
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


def _ToySettings(modes=3):
    return Kratos.Parameters("""{
        "input_dim"                  : %d,
        "context_dim"                : 1,
        "num_layers_decoder"         : 1,
        "num_heads"                  : 1,
        "dim_feedforward_scale"      : 2,
        "num_layers_context_encoder" : 1,
        "num_layers_input_encoder"   : 1,
        "num_layers_output_encoder"  : 1,
        "device"                     : "cpu"
    }""" % modes)


def _ToyTrajectories(samples=6, steps=8, modes=3, seed=0):
    """Damped rotating modes - a smooth, learnable q(t) family."""
    rng = numpy.random.default_rng(seed)
    trajectories = []
    for _ in range(samples):
        q0 = rng.uniform(0.5, 1.5, modes)
        t = numpy.arange(steps)[:, None]
        trajectories.append(q0 * numpy.exp(-0.1 * t) * numpy.cos(0.3 * t + q0))
    return numpy.stack(trajectories)  # (S, T, M)


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestRomTrajectoryDataset(KratosUnittest.TestCase):
    def test_ShapesAndDefaultContext(self):
        dataset = rom_temporal.CreateRomTrajectoryDataset(_ToyTrajectories())
        self.assertEqual(len(dataset), 6)
        z, context = dataset[2]
        self.assertEqual(tuple(z.shape), (8, 3))
        self.assertEqual(tuple(context.shape), (1, 1))
        self.assertEqual(float(context[0, 0]), 0.0)

    def test_ExplicitContexts(self):
        contexts = numpy.arange(6, dtype=float)[:, None]  # (S, 1)
        dataset = rom_temporal.CreateRomTrajectoryDataset(_ToyTrajectories(), contexts)
        _, context = dataset[4]
        self.assertEqual(tuple(context.shape), (1, 1))
        self.assertAlmostEqual(float(context[0, 0]), 4.0)

    def test_RaggedTrajectoriesRaise(self):
        with self.assertRaisesRegex(ValueError, "one \\(T, M\\) shape"):
            rom_temporal.CreateRomTrajectoryDataset(
                [numpy.zeros((8, 3)), numpy.zeros((5, 3))])

    def test_WrongContextShapeRaises(self):
        with self.assertRaisesRegex(ValueError, "contexts"):
            rom_temporal.CreateRomTrajectoryDataset(
                _ToyTrajectories(), numpy.zeros((2, 1)))


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestSequenceModelTraining(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_rom_temporal_model.pt")

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))
        if DistributedManager.is_initialized():
            DistributedManager.cleanup()

    def test_TwoEpochToyTraining(self):
        torch.manual_seed(0)
        model = rom_temporal.CreateSequenceModel(_ToySettings())
        dataset = rom_temporal.CreateRomTrajectoryDataset(_ToyTrajectories())
        history = rom_temporal.TrainRomTemporalModel(model, dataset, Kratos.Parameters("""{
            "epochs"     : 8,
            "batch_size" : 3,
            "device"     : "cpu",
            "seed"       : 0
        }"""))
        self.assertEqual(len(history), 8)
        self.assertTrue(all(numpy.isfinite(history)))
        self.assertLess(history[-1], history[0])
        self.assertFalse(model.training)

    def test_PredictRomTrajectoryShapes(self):
        model = rom_temporal.CreateSequenceModel(_ToySettings())
        single = rom_temporal.PredictRomTrajectory(model, numpy.zeros(3), steps=4)
        self.assertEqual(single.shape, (5, 3))
        prompt = rom_temporal.PredictRomTrajectory(model, numpy.zeros((2, 3)), steps=4)
        self.assertEqual(prompt.shape, (6, 3))
        self.assertTrue(numpy.isfinite(prompt).all())

    def test_SaveLoadRoundTrip(self):
        torch.manual_seed(1)
        settings = _ToySettings()
        model = rom_temporal.CreateSequenceModel(settings)
        rom_temporal.SaveRomTemporalModel(model, settings, self.checkpoint)

        restored, restored_settings = rom_temporal.LoadRomTemporalModel(self.checkpoint)
        self.assertEqual(restored_settings["input_dim"].GetInt(), 3)
        prompt = numpy.linspace(0.0, 1.0, 3)
        original = rom_temporal.PredictRomTrajectory(model.eval(), prompt, steps=3)
        reloaded = rom_temporal.PredictRomTrajectory(restored, prompt, steps=3)
        numpy.testing.assert_allclose(reloaded, original, atol=1e-6)

    def test_DistributedManagerAsDist(self):
        from KratosMultiphysics.PhysicsNeMoApplication.distributed import distributed_utils
        distributed_utils.InitializeDistributedManager(backend="gloo", port="29530")
        model = rom_temporal.CreateSequenceModel(_ToySettings())
        self.assertIsInstance(model.dist, DistributedManager)
        trajectory = rom_temporal.PredictRomTrajectory(model, numpy.zeros(3), steps=2)
        self.assertEqual(trajectory.shape, (3, 3))

    def test_InvalidInputDimRaises(self):
        with self.assertRaisesRegex(ValueError, "input_dim"):
            rom_temporal.CreateSequenceModel(Kratos.Parameters("""{ "device": "cpu" }"""))


if __name__ == '__main__':
    KratosUnittest.main()
