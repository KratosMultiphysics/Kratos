from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import physicsnemo.diffusion.samplers
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False

if have_physicsnemo:
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # physicsnemo.experimental import warning
            from physicsnemo.experimental.models import diffusion_unets  # noqa: F401
        have_unet3d = True
    except ImportError:
        have_unet3d = False
else:
    have_unet3d = False


def _WriteGrids(directory, steps, value_of_step, shape=(1, 8, 8, 2)):
    directory.mkdir(parents=True, exist_ok=True)
    for step in steps:
        numpy.savez(
            directory / f"grid_{step}.npz",
            grid=numpy.full(shape, value_of_step(step), dtype=numpy.float32),
            TIME=float(step), STEP=step,
            bounding_box=numpy.array([0.0, 0.0, -0.05, 1.0, 1.0, 0.05]))


def _TinyPreconditioner(seed=0):
    from physicsnemo.diffusion.preconditioners import EDMPrecondSuperResolution

    torch.manual_seed(seed)
    return EDMPrecondSuperResolution(
        img_resolution=8, img_in_channels=1, img_out_channels=1,
        model_type="SongUNet", model_channels=8, channel_mult=[1, 1],
        num_blocks=1, attn_resolutions=[])


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestGridPairDataset(KratosUnittest.TestCase):
    def setUp(self):
        self.lr_dir = Path("test_grid_pairs_lr")
        self.hr_dir = Path("test_grid_pairs_hr")

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.lr_dir))
        KratosUtilities.DeleteDirectoryIfExisting(str(self.hr_dir))

    def test_PairsMatchedByStep(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import CreateGridPairDataset

        _WriteGrids(self.lr_dir, (1, 2, 3), lambda s: 10.0 * s)
        _WriteGrids(self.hr_dir, (1, 2, 3), lambda s: 100.0 * s)
        dataset = CreateGridPairDataset(self.lr_dir, self.hr_dir, squeeze_axis=2)
        self.assertEqual(len(dataset), 3)
        condition, target = dataset[1]
        self.assertEqual(tuple(condition.shape), (1, 8, 8))
        self.assertTrue(bool((condition == 20.0).all()))
        self.assertTrue(bool((target == 200.0).all()))

    def test_StepMismatchRaises(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import CreateGridPairDataset

        _WriteGrids(self.lr_dir, (1, 2), lambda s: s)
        _WriteGrids(self.hr_dir, (1, 3), lambda s: s)
        with self.assertRaisesRegex(ValueError, "Step mismatch"):
            CreateGridPairDataset(self.lr_dir, self.hr_dir)


@KratosUnittest.skipUnless(have_torch and have_physicsnemo,
                           "Missing required python modules: torch, physicsnemo.")
class TestTrainDiffusionModel(KratosUnittest.TestCase):
    def test_LossIsFiniteOverEpochs(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils
        model = _TinyPreconditioner()
        torch.manual_seed(1)
        conditions = torch.randn(4, 1, 8, 8)
        targets = 2.0 * conditions + 0.1 * torch.randn(4, 1, 8, 8)
        dataset = torch.utils.data.TensorDataset(conditions, targets)

        history = diffusion_utils.TrainDiffusionModel(model, dataset, Kratos.Parameters("""{
            "epochs"     : 2,
            "batch_size" : 2,
            "device"     : "cpu",
            "seed"       : 3
        }"""))
        self.assertEqual(len(history), 2)
        self.assertTrue(all(numpy.isfinite(history)))
        self.assertFalse(model.training)  # left in eval mode

    def test_UnsupportedLossRaises(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils
        with self.assertRaisesRegex(ValueError, "diffusion loss"):
            diffusion_utils.TrainDiffusionModel(
                _TinyPreconditioner(), [], Kratos.Parameters("""{ "loss": "score_matching" }"""))

    def test_GenerateEnsembleShapes(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils
        model = _TinyPreconditioner().eval()
        ensemble = diffusion_utils.GenerateEnsemble(
            model, numpy.zeros((1, 8, 8)), Kratos.Parameters("""{
                "num_samples" : 2,
                "num_steps"   : 4,
                "seed"        : 5
            }"""))
        self.assertEqual(ensemble.shape, (2, 1, 8, 8))
        self.assertTrue(numpy.isfinite(ensemble).all())
        # two independent samples of a diffusion model differ
        self.assertGreater(numpy.abs(ensemble[0] - ensemble[1]).max(), 0.0)


@KratosUnittest.skipUnless(have_torch and have_physicsnemo,
                           "Missing required python modules: torch, physicsnemo.")
class TestDiffusionInferenceProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_diffusion_model.mdlus")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE, Kratos.NODAL_ERROR))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + node.Y)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def test_EnsembleMeanAndUncertaintyAreWritten(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import diffusion_inference_process
        _TinyPreconditioner().save(str(self.checkpoint))
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"    : "Main",
                "model_settings"     : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "input_fields"       : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"      : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "uncertainty_fields" : [ { "variable_name" : "NODAL_ERROR", "data_location" : "node_historical" } ],
                "grid_shape"         : [8, 8, 2],
                "squeeze_axis"       : 2,
                "sampler_settings"   : { "num_samples" : 2, "num_steps" : 4, "seed" : 7 }
            }
        }""" % self.checkpoint)
        process = diffusion_inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        mean = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        spread = numpy.array([
            node.GetSolutionStepValue(Kratos.NODAL_ERROR) for node in self.model_part.Nodes])
        self.assertTrue(numpy.isfinite(mean).all())
        self.assertGreater(numpy.abs(mean).max(), 0.0)
        self.assertTrue(numpy.isfinite(spread).all())
        self.assertTrue((spread >= 0.0).all())
        self.assertGreater(spread.max(), 0.0)

    def test_OutputNormalizationFromTheCardIsApplied(self):
        """Same seeded sampler with and without a card: the mean is scaled
        AND shifted, the spread scaled only - a spread shifted by the
        training mean would be the shared-hook mistake."""
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import diffusion_inference_process
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
        _TinyPreconditioner().save(str(self.checkpoint))

        def Run():
            settings = Kratos.Parameters("""{
                "Parameters": {
                    "model_part_name"    : "Main",
                    "model_settings"     : {
                        "checkpoint_file" : "%s",
                        "checkpoint_type" : "physicsnemo",
                        "device"          : "cpu"
                    },
                    "input_fields"       : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                    "output_fields"      : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                    "uncertainty_fields" : [ { "variable_name" : "NODAL_ERROR", "data_location" : "node_historical" } ],
                    "grid_shape"         : [8, 8, 2],
                    "squeeze_axis"       : 2,
                    "sampler_settings"   : { "num_samples" : 2, "num_steps" : 4, "seed" : 7 }
                }
            }""" % self.checkpoint)
            process = diffusion_inference_process.Factory(settings, self.model)
            self.model_part.ProcessInfo[Kratos.STEP] = 1
            process.ExecuteFinalizeSolutionStep()
            return (numpy.array([node.GetSolutionStepValue(Kratos.TEMPERATURE)
                                 for node in self.model_part.Nodes]),
                    numpy.array([node.GetSolutionStepValue(Kratos.NODAL_ERROR)
                                 for node in self.model_part.Nodes]))

        raw_mean, raw_spread = Run()

        # input side first: the card standardizing the raw condition must
        # reproduce the pre-standardized-field run exactly (same seed)
        in_mean, in_std = 0.5, 2.0
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, (node.X + node.Y - in_mean) / in_std)
        reference_mean, reference_spread = Run()
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + node.Y)
        model_registry.SaveModelCard(self.checkpoint, {
            "input_normalization": {"type": "mean_std", "mean": [in_mean], "std": [in_std]}})
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")
        card_mean, card_spread = Run()
        numpy.testing.assert_allclose(card_mean, reference_mean, rtol=1e-6, atol=1e-9)
        numpy.testing.assert_allclose(card_spread, reference_spread, rtol=1e-6, atol=1e-9)
        self.assertGreater(numpy.abs(reference_mean - raw_mean).max(), 0.0)

        mean, std = 50.0, 4.0
        model_registry.SaveModelCard(self.checkpoint, {
            "output_normalization": {"type": "mean_std", "mean": [mean], "std": [std]}})
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")
        card_mean, card_spread = Run()
        numpy.testing.assert_allclose(card_mean, std * raw_mean + mean, rtol=1e-6, atol=1e-9)
        numpy.testing.assert_allclose(card_spread, std * raw_spread, rtol=1e-6, atol=1e-9)

    def test_InvalidSqueezeAxisRaises(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import diffusion_inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "squeeze_axis"    : 5
            }
        }""")
        with self.assertRaisesRegex(ValueError, "squeeze_axis"):
            diffusion_inference_process.Factory(settings, self.model)


@KratosUnittest.skipUnless(have_torch and have_physicsnemo,
                           "Missing required python modules: torch, physicsnemo.")
class TestDitDenoiser(KratosUnittest.TestCase):
    """physicsnemo.models.dit.DiT as the denoiser, via WrapDenoiser and the
    process's "denoiser_interface": "dit"."""

    def setUp(self):
        self.checkpoint = Path("test_diffusion_dit.mdlus")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE, Kratos.NODAL_ERROR))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + node.Y)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    @staticmethod
    def _TinyDit(seed=0):
        from physicsnemo.models.dit import DiT

        torch.manual_seed(seed)
        # in_channels = 1 latent + 1 condition (concatenated by the wrapper)
        return DiT(input_size=8, in_channels=2, out_channels=1, patch_size=4,
                   hidden_size=32, depth=1, num_heads=2)

    def test_WrapDenoiserContract(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils
        wrapper = diffusion_utils.WrapDenoiser(self._TinyDit(), "dit")
        self.assertEqual(wrapper.img_out_channels, 1)

        x = torch.randn(2, 1, 8, 8)
        img_lr = torch.randn(2, 1, 8, 8)
        with torch.no_grad():
            out_scalar_sigma = wrapper(x, img_lr, torch.tensor(0.5))
            out_batch_sigma = wrapper(x, img_lr, torch.tensor([0.5, 1.0]))
        self.assertEqual(list(out_scalar_sigma.shape), [2, 1, 8, 8])
        self.assertEqual(list(out_batch_sigma.shape), [2, 1, 8, 8])

        with self.assertRaisesRegex(ValueError, "denoiser interface"):
            diffusion_utils.WrapDenoiser(self._TinyDit(), "songunet")

    def test_DitThroughDiffusionProcess(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import diffusion_inference_process
        self._TinyDit().save(str(self.checkpoint))
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"    : "Main",
                "model_settings"     : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "denoiser_interface" : "dit",
                "input_fields"       : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"      : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "uncertainty_fields" : [ { "variable_name" : "NODAL_ERROR", "data_location" : "node_historical" } ],
                "grid_shape"         : [8, 8, 2],
                "squeeze_axis"       : 2,
                "sampler_settings"   : { "num_samples" : 2, "num_steps" : 4, "seed" : 7, "output_channels" : 1 }
            }
        }""" % self.checkpoint)
        process = diffusion_inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        means = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        spreads = numpy.array([
            node.GetSolutionStepValue(Kratos.NODAL_ERROR) for node in self.model_part.Nodes])
        self.assertTrue(numpy.isfinite(means).all())
        self.assertTrue(numpy.isfinite(spreads).all())
        self.assertTrue((spreads >= 0.0).all())
        self.assertGreater(spreads.max(), 0.0)

    def test_UnknownDenoiserInterfaceRaises(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import diffusion_inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"    : "Main",
                "denoiser_interface" : "songunet",
                "input_fields"       : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"      : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""")
        with self.assertRaisesRegex(ValueError, "denoiser interface"):
            diffusion_inference_process.Factory(settings, self.model)


@KratosUnittest.skipUnless(have_torch and have_unet3d,
                           "Missing required python modules: torch, physicsnemo (with "
                           "experimental DiffusionUNet3D).")
class TestUnet3dDenoiser(KratosUnittest.TestCase):
    """The volumetric physicsnemo.experimental.models.diffusion_unets.DiffusionUNet3D
    as the denoiser, via WrapDenoiser("unet3d") and the process's
    "denoiser_interface": "unet3d" - full (C, D, H, W) grids, no thin-axis squeeze."""

    def setUp(self):
        self.checkpoint = Path("test_diffusion_unet3d.mdlus")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE, Kratos.NODAL_ERROR))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + node.Y + node.Z)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    @staticmethod
    def _TinyUnet3d(seed=0):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # physicsnemo.experimental import warning
            from physicsnemo.experimental.models.diffusion_unets import DiffusionUNet3D

        torch.manual_seed(seed)
        # x_channels = the 1 predicted channel; the condition is NOT
        # concatenated - it feeds vol_cond_channels natively
        return DiffusionUNet3D(
            x_channels=1, vol_cond_channels=1, num_levels=2, model_channels=8,
            channel_mult=[1, 1], num_blocks=1, attention_levels=[],
            bottleneck_attention=False, dropout=0.0)

    def test_WrapDenoiserContract(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils
        wrapper = diffusion_utils.WrapDenoiser(self._TinyUnet3d(), "unet3d")
        self.assertEqual(wrapper.img_out_channels, 1)

        x = torch.randn(2, 1, 4, 4, 4)
        img_lr = torch.randn(2, 1, 4, 4, 4)
        with torch.no_grad():
            out_scalar_sigma = wrapper(x, img_lr, torch.tensor(0.5))
            out_batch_sigma = wrapper(x, img_lr, torch.tensor([0.5, 1.0]))
        self.assertEqual(list(out_scalar_sigma.shape), [2, 1, 4, 4, 4])
        self.assertEqual(list(out_batch_sigma.shape), [2, 1, 4, 4, 4])

        # a squeezed planar latent cannot feed the volumetric interface
        with self.assertRaisesRegex(ValueError, "volumetric"):
            wrapper(torch.randn(2, 1, 8, 8), torch.randn(2, 1, 8, 8), torch.tensor(0.5))
        # a model without x_channels needs out_channels set explicitly
        with self.assertRaisesRegex(ValueError, "x_channels"):
            diffusion_utils.WrapDenoiser(torch.nn.Conv3d(1, 1, 1), "unet3d")

    def test_VolumetricTraining(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils
        wrapper = diffusion_utils.WrapDenoiser(self._TinyUnet3d(), "unet3d")
        torch.manual_seed(1)
        conditions = torch.randn(4, 1, 4, 4, 4)
        targets = 2.0 * conditions + 0.1 * torch.randn(4, 1, 4, 4, 4)
        dataset = torch.utils.data.TensorDataset(conditions, targets)

        history = diffusion_utils.TrainDiffusionModel(wrapper, dataset, Kratos.Parameters("""{
            "epochs"     : 2,
            "batch_size" : 2,
            "device"     : "cpu",
            "seed"       : 3
        }"""))
        self.assertEqual(len(history), 2)
        self.assertTrue(all(numpy.isfinite(history)))
        self.assertFalse(wrapper.training)  # left in eval mode

        # the CorrDiff losses are 2D-only upstream and must say so
        with self.assertRaisesRegex(ValueError, "2D-only"):
            diffusion_utils.TrainDiffusionModel(wrapper, dataset, Kratos.Parameters("""{
                "loss" : "regression"
            }"""))

    def test_Unet3dThroughDiffusionProcess(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import diffusion_inference_process
        self._TinyUnet3d().save(str(self.checkpoint))
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"    : "Main",
                "model_settings"     : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "denoiser_interface" : "unet3d",
                "input_fields"       : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"      : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "uncertainty_fields" : [ { "variable_name" : "NODAL_ERROR", "data_location" : "node_historical" } ],
                "grid_shape"         : [4, 4, 4],
                "sampler_settings"   : { "num_samples" : 2, "num_steps" : 4, "seed" : 7, "output_channels" : 1 }
            }
        }""" % self.checkpoint)
        process = diffusion_inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        means = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        spreads = numpy.array([
            node.GetSolutionStepValue(Kratos.NODAL_ERROR) for node in self.model_part.Nodes])
        self.assertTrue(numpy.isfinite(means).all())
        self.assertGreater(numpy.abs(means).max(), 0.0)
        self.assertTrue(numpy.isfinite(spreads).all())
        self.assertTrue((spreads >= 0.0).all())
        self.assertGreater(spreads.max(), 0.0)

    def test_SqueezeAxisRejectedForUnet3d(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import diffusion_inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"    : "Main",
                "denoiser_interface" : "unet3d",
                "squeeze_axis"       : 2,
                "input_fields"       : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"      : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""")
        with self.assertRaisesRegex(ValueError, "squeeze_axis"):
            diffusion_inference_process.Factory(settings, self.model)


if __name__ == '__main__':
    KratosUnittest.main()
