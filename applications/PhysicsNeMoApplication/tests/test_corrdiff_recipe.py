"""Tests for the CorrDiff two-stage recipe (regression + residual
diffusion) and the FWI-style inversion recipe on the diffusion bridge."""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils
try:
    import torch
    from physicsnemo.diffusion.preconditioners import EDMPrecondSuperResolution
    from physicsnemo.models.diffusion_unets import CorrDiffRegressionUNet
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


def _TinyRegression(in_channels=2, out_channels=1, seed=0):
    torch.manual_seed(seed)
    return CorrDiffRegressionUNet(
        img_resolution=8, img_in_channels=in_channels, img_out_channels=out_channels,
        model_type="SongUNet", model_channels=8, channel_mult=[1, 1],
        num_blocks=1, attn_resolutions=[])


def _TinyDenoiser(in_channels=2, out_channels=1, seed=1):
    torch.manual_seed(seed)
    return EDMPrecondSuperResolution(
        img_resolution=8, img_in_channels=in_channels, img_out_channels=out_channels,
        model_type="SongUNet", model_channels=8, channel_mult=[1, 1],
        num_blocks=1, attn_resolutions=[])


def _TinyResidualDenoiser(in_channels=2, out_channels=1, seed=1):
    """ResidualLoss passes embedding_selector/global_index kwargs, so the
    residual-stage denoiser must wrap SongUNetPosEmbd - whose
    N_grid_channels (default 4) positional channels count toward
    img_in_channels, exactly as upstream CorrDiff configs size them."""
    torch.manual_seed(seed)
    return EDMPrecondSuperResolution(
        img_resolution=8, img_in_channels=in_channels + 4, img_out_channels=out_channels,
        model_type="SongUNetPosEmbd", model_channels=8, channel_mult=[1, 1],
        num_blocks=1, attn_resolutions=[])


def _WindLikeDataset(samples=8, seed=0):
    """Synthetic downscaling pairs: condition = (smooth field, fixed
    topography), target = a deterministic map of both + small noise."""
    rng = numpy.random.default_rng(seed)
    yy, xx = numpy.meshgrid(numpy.linspace(0, 1, 8), numpy.linspace(0, 1, 8), indexing="ij")
    topography = numpy.sin(2 * numpy.pi * xx) * numpy.cos(numpy.pi * yy)  # static channel

    conditions, targets = [], []
    for _ in range(samples):
        smooth = rng.standard_normal((8, 8))
        for _ in range(3):
            smooth = (smooth + numpy.roll(smooth, 1, 0) + numpy.roll(smooth, -1, 0)
                      + numpy.roll(smooth, 1, 1) + numpy.roll(smooth, -1, 1)) / 5.0
        conditions.append(numpy.stack([smooth, topography]))
        targets.append((0.7 * smooth + 0.3 * topography
                        + 0.02 * rng.standard_normal((8, 8)))[None])
    return (torch.tensor(numpy.stack(conditions), dtype=torch.float32),
            torch.tensor(numpy.stack(targets), dtype=torch.float32))


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python modules: torch, physicsnemo.")
class TestCorrDiffTwoStage(KratosUnittest.TestCase):
    def test_RegressionStageLearnsTheMean(self):
        conditions, targets = _WindLikeDataset()
        dataset = torch.utils.data.TensorDataset(conditions, targets)
        regression = _TinyRegression()

        history = diffusion_utils.TrainDiffusionModel(regression, dataset, Kratos.Parameters("""{
            "epochs"        : 30,
            "batch_size"    : 4,
            "learning_rate" : 2e-3,
            "loss"          : "regression",
            "device"        : "cpu",
            "seed"          : 0
        }"""))
        self.assertLess(history[-1], history[0])

        # the learned mean must beat the zero baseline on the training set
        errors, baselines = [], []
        for i in range(len(conditions)):
            mean = diffusion_utils.RunRegressionMean(regression, conditions[i].numpy())
            errors.append(float(numpy.mean((mean - targets[i].numpy()) ** 2)))
            baselines.append(float(numpy.mean(targets[i].numpy() ** 2)))
        self.assertLess(numpy.mean(errors), 0.5 * numpy.mean(baselines))

    def test_ResidualStageRunsWithFrozenRegression(self):
        conditions, targets = _WindLikeDataset()
        dataset = torch.utils.data.TensorDataset(conditions, targets)
        regression = _TinyRegression()
        denoiser = _TinyResidualDenoiser()

        regression_history, diffusion_history = diffusion_utils.TrainCorrDiffPair(
            regression, denoiser, dataset, Kratos.Parameters("""{
                "epochs"            : 3,
                "regression_epochs" : 10,
                "batch_size"        : 4,
                "learning_rate"     : 1e-3,
                "device"            : "cpu",
                "seed"              : 0
            }"""))
        self.assertEqual(len(regression_history), 10)
        self.assertEqual(len(diffusion_history), 3)
        self.assertTrue(all(numpy.isfinite(diffusion_history)))
        # the regression stage stayed frozen during stage 2
        self.assertFalse(any(p.requires_grad for p in regression.parameters()))

    def test_ResidualWithoutRegressionRaises(self):
        conditions, targets = _WindLikeDataset(samples=2)
        dataset = torch.utils.data.TensorDataset(conditions, targets)
        with self.assertRaisesRegex(ValueError, "regression_model"):
            diffusion_utils.TrainDiffusionModel(
                _TinyDenoiser(), dataset, Kratos.Parameters('{"loss": "residual", "epochs": 1}'))

    def test_UnknownLossRaises(self):
        conditions, targets = _WindLikeDataset(samples=2)
        dataset = torch.utils.data.TensorDataset(conditions, targets)
        with self.assertRaisesRegex(ValueError, "diffusion loss"):
            diffusion_utils.TrainDiffusionModel(
                _TinyDenoiser(), dataset, Kratos.Parameters('{"loss": "score_matching"}'))

    def test_RunRegressionMeanValidation(self):
        class NoChannels(torch.nn.Module):
            def forward(self, x, img_lr):
                return x

        with self.assertRaisesRegex(ValueError, "img_out_channels"):
            diffusion_utils.RunRegressionMean(NoChannels(), numpy.zeros((2, 8, 8)))


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python modules: torch, physicsnemo.")
class TestRegressionSettingsThroughProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.regression_checkpoint = Path("test_corrdiff_regression.mdlus")
        self.denoiser_checkpoint = Path("test_corrdiff_denoiser.mdlus")
        from test_grid_bridge import CreateStructuredTetModelPart
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE, Kratos.NODAL_PAUX))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + node.Y)

    def tearDown(self):
        for path in (self.regression_checkpoint, self.denoiser_checkpoint):
            KratosUtilities.DeleteFileIfExisting(str(path))
            KratosUtilities.DeleteFileIfExisting(str(path) + ".card.json")

    def _CreateProcess(self, output_variable, with_regression):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import diffusion_inference_process
        regression_block = ("""
                "regression_settings" : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },""" % self.regression_checkpoint) if with_regression else ""
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"    : "Main",
                "model_settings"     : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },%s
                "input_fields"       : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
                "output_fields"      : [ { "variable_name" : "%s", "data_location" : "node_historical" } ],
                "grid_shape"         : [8, 8, 2],
                "squeeze_axis"       : 2,
                "sampler_settings"   : { "num_samples" : 2, "num_steps" : 4, "seed" : 7 }
            }
        }""" % (self.denoiser_checkpoint, regression_block, output_variable))
        return diffusion_inference_process.Factory(settings, self.model)

    def test_RegressionMeanShiftsTheEnsemble(self):
        _TinyRegression(in_channels=1).save(str(self.regression_checkpoint))
        _TinyDenoiser(in_channels=1).save(str(self.denoiser_checkpoint))

        # same denoiser + same sampler seed: the with-regression run differs
        # from the plain run by EXACTLY the scattered regression mean
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        self._CreateProcess("TEMPERATURE", with_regression=False).ExecuteFinalizeSolutionStep()
        self._CreateProcess("NODAL_PAUX", with_regression=True).ExecuteFinalizeSolutionStep()

        from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
        # recompute the regression mean on the process's exact condition grid
        bounding_box = grid_bridge.ComputeBoundingBox(self.model_part)
        condition, _ = grid_bridge.SampleFieldsOnGrid(
            self.model_part, [("PRESSURE", "node_historical")], (8, 8, 2), bounding_box)
        condition = condition.mean(axis=3)  # squeeze_axis 2 -> mean over axis 1+2
        regression, _ = model_registry.LoadModel(Kratos.Parameters("""{
            "checkpoint_file" : "%s",
            "checkpoint_type" : "physicsnemo",
            "device"          : "cpu"
        }""" % self.regression_checkpoint))
        mean_grid = diffusion_utils.RunRegressionMean(regression, condition)
        # duplicate across the squeezed thin axis and scatter like the process
        mean_grid = numpy.repeat(numpy.expand_dims(mean_grid, 3), 2, axis=3)
        for node in self.model_part.Nodes:
            node.SetValue(Kratos.NODAL_ERROR, 0.0)
        grid_bridge.ScatterGridToNodes(
            mean_grid, bounding_box, self.model_part,
            [("NODAL_ERROR", "node_non_historical")])

        for node in self.model_part.Nodes:
            difference = (node.GetSolutionStepValue(Kratos.NODAL_PAUX)
                          - node.GetSolutionStepValue(Kratos.TEMPERATURE))
            self.assertAlmostEqual(difference, node.GetValue(Kratos.NODAL_ERROR), places=5)


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python modules: torch, physicsnemo.")
class TestFwiInversionRecipe(KratosUnittest.TestCase):
    """Inversion by conditional diffusion: layered-earth property grids
    conditioned on sparse borehole observations + an observation mask -
    nothing beyond the shipped diffusion bridge is needed."""

    def test_LayeredEarthInversion(self):
        rng = numpy.random.default_rng(0)
        conditions, targets = [], []
        for _ in range(8):
            # layered earth: piecewise-constant-by-row property
            depths = numpy.sort(rng.integers(1, 7, size=2))
            values = rng.uniform(0.2, 1.0, size=3)
            profile = numpy.empty(8)
            profile[:depths[0]] = values[0]
            profile[depths[0]:depths[1]] = values[1]
            profile[depths[1]:] = values[2]
            property_grid = numpy.tile(profile[:, None], (1, 8))
            # observations: two random borehole columns + the mask channel
            observed = numpy.zeros((8, 8))
            mask = numpy.zeros((8, 8))
            for column in rng.choice(8, size=2, replace=False):
                observed[:, column] = property_grid[:, column]
                mask[:, column] = 1.0
            conditions.append(numpy.stack([observed, mask]))
            targets.append(property_grid[None])

        dataset = torch.utils.data.TensorDataset(
            torch.tensor(numpy.stack(conditions), dtype=torch.float32),
            torch.tensor(numpy.stack(targets), dtype=torch.float32))

        denoiser = _TinyDenoiser(in_channels=2)
        history = diffusion_utils.TrainDiffusionModel(denoiser, dataset, Kratos.Parameters("""{
            "epochs"        : 8,
            "batch_size"    : 4,
            "learning_rate" : 1e-3,
            "device"        : "cpu",
            "seed"          : 0
        }"""))
        # the EDM loss samples sigma stochastically - per-epoch values are
        # noisy, so assert on the best epoch, not the last
        self.assertTrue(numpy.isfinite(history).all())
        self.assertLess(min(history[1:]), history[0])

        ensemble = diffusion_utils.GenerateEnsemble(
            denoiser, conditions[0], Kratos.Parameters("""{
                "num_samples" : 4, "num_steps" : 8, "seed" : 0
            }"""))
        self.assertEqual(ensemble.shape, (4, 1, 8, 8))
        self.assertTrue(numpy.isfinite(ensemble).all())
        # the ensemble spread is the inversion uncertainty - it must exist
        self.assertGreater(ensemble.std(axis=0).max(), 0.0)


if __name__ == '__main__':
    KratosUnittest.main()
