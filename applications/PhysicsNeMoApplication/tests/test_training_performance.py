"""Tests for TrainModel's performance layer: AMP, static capture, gradient
clipping, the Muon optimizer, learning-rate schedules, profiling, the
LaunchLogger, and resumable checkpoints."""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.models.mlp import FullyConnected
    from physicsnemo.utils import StaticCaptureTraining  # noqa: F401
    from physicsnemo.optim import Muon  # noqa: F401
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


def _Model(seed=0):
    torch.manual_seed(seed)
    return FullyConnected(in_features=3, out_features=1, num_layers=2, layer_size=16)


def _Dataset(samples=64, seed=0):
    generator = torch.Generator().manual_seed(seed)
    inputs = torch.randn(samples, 3, generator=generator)
    return torch.utils.data.TensorDataset(inputs, inputs.sum(dim=1, keepdim=True))


def _Settings(epochs=6, **performance):
    import json

    settings = Kratos.Parameters("""{
        "epochs"        : %d,
        "batch_size"    : 16,
        "learning_rate" : 1e-2,
        "device"        : "cpu",
        "shuffle"       : false,
        "seed"          : 0
    }""" % epochs)
    if performance:
        settings.AddValue("performance", Kratos.Parameters(json.dumps(performance)))
    return settings


@KratosUnittest.skipUnless(have_torch and have_physicsnemo,
                           "Missing required python modules: torch, physicsnemo.")
class TestTrainingPerformanceLayer(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoints = Path("test_training_checkpoints")
        self.profile = Path("test_training_profile.json")

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting(str(self.checkpoints))
        kratos_utils.DeleteFileIfExisting(str(self.profile))

    def test_AnEmptyBlockChangesNothing(self):
        plain = training_utils.TrainModel(_Model(), _Dataset(), _Settings())
        empty = training_utils.TrainModel(_Model(), _Dataset(), _Settings(**{}))
        numpy.testing.assert_allclose(plain, empty, rtol=0, atol=0)

    def test_AmpInBfloat16TrainsOnCpu(self):
        history = training_utils.TrainModel(
            _Model(), _Dataset(), _Settings(epochs=10, amp=True, amp_dtype="bfloat16"))
        self.assertTrue(numpy.isfinite(history).all())
        self.assertLess(history[-1], history[0])

    def test_StaticCaptureTrainsAPhysicsNemoModule(self):
        history = training_utils.TrainModel(
            _Model(), _Dataset(), _Settings(epochs=10, static_capture=True))
        self.assertTrue(numpy.isfinite(history).all())
        self.assertLess(history[-1], history[0])

    def test_StaticCaptureRefusesAPlainTorchModule(self):
        """physicsnemo's StaticCaptureTraining accepts only its own Module
        type; the message says how to get one."""
        with self.assertRaisesRegex(ValueError, "physicsnemo Module"):
            training_utils.TrainModel(
                torch.nn.Linear(3, 1), _Dataset(), _Settings(epochs=1, static_capture=True))

    def test_GradientClippingRuns(self):
        history = training_utils.TrainModel(
            _Model(), _Dataset(), _Settings(epochs=4, gradient_clip_norm=0.5))
        self.assertTrue(numpy.isfinite(history).all())

    def test_MuonOnTheMatricesAndAdamOnTheRest(self):
        settings = _Settings(epochs=10)
        settings.AddEmptyValue("optimizer").SetString("muon")
        history = training_utils.TrainModel(_Model(), _Dataset(), settings)
        self.assertLess(history[-1], history[0])

    def test_MuonNeedsAMatrixToOptimize(self):
        """physicsnemo's Muon rejects parameters of rank < 2 outright, so a
        model with none has nothing for it to do."""
        class BiasOnly(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.bias = torch.nn.Parameter(torch.zeros(1))

            def forward(self, x):
                return x[..., :1] * 0.0 + self.bias

        settings = _Settings(epochs=1)
        settings.AddEmptyValue("optimizer").SetString("muon")
        with self.assertRaisesRegex(ValueError, "rank >= 2"):
            training_utils.TrainModel(BiasOnly(), _Dataset(), settings)

    def test_ProfileWritesAChromeTrace(self):
        training_utils.TrainModel(
            _Model(), _Dataset(samples=16),
            _Settings(epochs=1, profile=True, profile_output=str(self.profile)))
        self.assertTrue(self.profile.exists())
        self.assertGreater(self.profile.stat().st_size, 0)

    def test_LaunchLoggerOnTheConsole(self):
        history = training_utils.TrainModel(
            _Model(), _Dataset(samples=16), _Settings(epochs=2, launch_logger=True))
        self.assertEqual(len(history), 2)

    def test_ResumingReproducesTheUninterruptedRun(self):
        """The property that makes a checkpoint resumable rather than merely
        loadable: model, optimizer AND scheduler state all come back, so the
        second half of a split run is the second half of the whole one."""
        schedule = {"type": "step", "step_size": 1, "gamma": 0.7}
        whole = training_utils.TrainModel(
            _Model(), _Dataset(), _Settings(epochs=4, scheduler=schedule))

        first_half = training_utils.TrainModel(
            _Model(), _Dataset(), _Settings(
                epochs=2, scheduler=schedule, checkpoint_directory=str(self.checkpoints),
                checkpoint_interval=1))
        second_half = training_utils.TrainModel(
            _Model(), _Dataset(), _Settings(
                epochs=4, scheduler=schedule, checkpoint_directory=str(self.checkpoints),
                checkpoint_interval=1, resume=True))

        self.assertEqual(len(first_half), 2)
        self.assertEqual(len(second_half), 2)   # only the epochs run in THIS call
        numpy.testing.assert_allclose(first_half, whole[:2], rtol=1e-10)
        numpy.testing.assert_allclose(second_half, whole[2:], rtol=1e-6)

    def test_CheckpointFilesAreWrittenPerInterval(self):
        training_utils.TrainModel(
            _Model(), _Dataset(samples=16), _Settings(
                epochs=3, checkpoint_directory=str(self.checkpoints), checkpoint_interval=2))
        epochs_saved = sorted(
            int(path.name.split(".")[-2]) for path in self.checkpoints.glob("checkpoint.*.pt"))
        self.assertEqual(epochs_saved, [2, 3])  # every 2 epochs, and the last one

    def test_ResumingFromAnEmptyDirectoryStartsFresh(self):
        history = training_utils.TrainModel(
            _Model(), _Dataset(samples=16), _Settings(
                epochs=2, checkpoint_directory=str(self.checkpoints), resume=True))
        self.assertEqual(len(history), 2)

    def test_InvalidPerformanceSettingsRaise(self):
        for block, message in ((dict(resume=True), "checkpoint_directory"),
                               (dict(amp_dtype="float8"), "amp_dtype"),
                               (dict(logger_backend="tensorboard"), "logger_backend"),
                               (dict(scheduler={"type": "polynomial"}), "scheduler")):
            with self.assertRaisesRegex(ValueError, message):
                training_utils.TrainModel(_Model(), _Dataset(samples=8), _Settings(epochs=1, **block))


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestTrainingSchedulers(KratosUnittest.TestCase):
    @staticmethod
    def _Optimizer():
        return torch.optim.SGD([torch.nn.Parameter(torch.zeros(2, 2))], lr=1.0)

    def _Trajectory(self, scheduler_settings, epochs):
        optimizer = self._Optimizer()
        scheduler = training_utils._BuildScheduler(optimizer, scheduler_settings, epochs)
        rates = []
        for _ in range(epochs):
            rates.append(optimizer.param_groups[0]["lr"])
            optimizer.step()
            scheduler.step()
        return rates

    def test_StepHalvesOnSchedule(self):
        rates = self._Trajectory({"type": "step", "step_size": 2, "gamma": 0.5, "t_max": 0}, 6)
        numpy.testing.assert_allclose(rates, [1.0, 1.0, 0.5, 0.5, 0.25, 0.25])

    def test_CosineAnnealsToZeroOverTheRun(self):
        rates = self._Trajectory({"type": "cosine", "step_size": 1, "gamma": 0.1, "t_max": 0}, 4)
        self.assertAlmostEqual(rates[0], 1.0)
        self.assertLess(rates[-1], rates[0])
        optimizer = self._Optimizer()
        scheduler = training_utils._BuildScheduler(
            optimizer, {"type": "cosine", "step_size": 1, "gamma": 0.1, "t_max": 0}, 4)
        for _ in range(4):
            optimizer.step()
            scheduler.step()
        self.assertAlmostEqual(optimizer.param_groups[0]["lr"], 0.0, places=12)

    def test_NoneIsNoScheduler(self):
        self.assertIsNone(training_utils._BuildScheduler(
            self._Optimizer(), {"type": "none", "step_size": 1, "gamma": 0.1, "t_max": 0}, 3))


if __name__ == '__main__':
    KratosUnittest.main()
