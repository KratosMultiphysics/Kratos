import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import superresolution_process
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
sys.path.insert(0, str(Path(__file__).parent))
from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.models.srrn import SRResNet
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestSuperResolutionProcess(KratosUnittest.TestCase):
    """Uses a scripted trilinear upsampler as the 'model': exact on linear
    fields, so the fine-mesh result is analytically checkable."""

    def setUp(self):
        self.checkpoint = Path("test_superresolution_upsampler.pt")

        class Upsampler(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.upsample = torch.nn.Upsample(scale_factor=2, mode="trilinear", align_corners=True)

            def forward(self, x):
                return self.upsample(x)

        torch.jit.script(Upsampler()).save(str(self.checkpoint))

        self.model = Kratos.Model()
        self.coarse = CreateStructuredTetModelPart(self.model, "Coarse", divisions=3)
        self.fine = CreateStructuredTetModelPart(self.model, "Fine", divisions=5,
                                                 historical_variables=(Kratos.TEMPERATURE,))
        for node in self.coarse.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0 + 2.0 * node.X + 3.0 * node.Y - node.Z)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def _CreateProcess(self, output_interval=1):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "coarse_model_part_name" : "Coarse",
                "fine_model_part_name"   : "Fine",
                "model_settings"         : {
                    "checkpoint_file" : "test_superresolution_upsampler.pt",
                    "device"          : "cpu"
                },
                "input_fields"           : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"          : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "coarse_grid_shape"      : [6, 6, 6],
                "output_interval"        : 1
            }
        }""")
        settings["Parameters"]["output_interval"].SetInt(output_interval)
        return superresolution_process.Factory(settings, self.model)

    def test_LinearFieldSuperresolvedExactly(self):
        process = self._CreateProcess()
        self.coarse.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        for node in self.fine.Nodes:
            expected = 1.0 + 2.0 * node.X + 3.0 * node.Y - node.Z
            self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.TEMPERATURE), expected, places=8)

    def _CreateTwoChannelProcess(self, card):
        """Two channels with different scalings: the card's per-channel
        vectors must run along axis 0 (the grid layout), so applying them
        along the last, spatial axis raises instead of passing by
        broadcast. PRESSURE = 1 + 2x, TEMPERATURE = 3y - z on the coarse
        part; the upsampler is exact on both."""
        coarse = CreateStructuredTetModelPart(
            self.model, "Coarse2", divisions=3,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        fine = CreateStructuredTetModelPart(
            self.model, "Fine2", divisions=5,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        for node in coarse.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0 + 2.0 * node.X)
            node.SetSolutionStepValue(Kratos.TEMPERATURE, 3.0 * node.Y - node.Z)
        model_registry.SaveModelCard(self.checkpoint, card)
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")

        settings = Kratos.Parameters("""{
            "Parameters": {
                "coarse_model_part_name" : "Coarse2",
                "fine_model_part_name"   : "Fine2",
                "model_settings"         : {
                    "checkpoint_file" : "test_superresolution_upsampler.pt",
                    "device"          : "cpu"
                },
                "input_fields"           : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" },
                                             { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "output_fields"          : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" },
                                             { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "coarse_grid_shape"      : [6, 6, 6]
            }
        }""")
        process = superresolution_process.Factory(settings, self.model)
        coarse.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        return fine

    def test_OutputNormalizationFromTheCardIsApplied(self):
        mean, std = [-4.0, 100.0], [2.5, 0.1]
        fine = self._CreateTwoChannelProcess({
            "output_normalization": {"type": "mean_std", "mean": mean, "std": std}})
        for node in fine.Nodes:
            self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.TEMPERATURE),
                                   std[0] * (1.0 + 2.0 * node.X) + mean[0], places=8)
            self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.PRESSURE),
                                   std[1] * (3.0 * node.Y - node.Z) + mean[1], places=8)

    def test_InputNormalizationFromTheCardIsApplied(self):
        mean, std = [1.0, -2.0], [2.0, 0.5]
        fine = self._CreateTwoChannelProcess({
            "input_normalization": {"type": "mean_std", "mean": mean, "std": std}})
        for node in fine.Nodes:
            self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.TEMPERATURE),
                                   ((1.0 + 2.0 * node.X) - mean[0]) / std[0], places=8)
            self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.PRESSURE),
                                   ((3.0 * node.Y - node.Z) - mean[1]) / std[1], places=8)

    def test_IntervalGating(self):
        process = self._CreateProcess(output_interval=2)
        self.coarse.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        self.assertEqual(self.fine.GetNode(1).GetSolutionStepValue(Kratos.TEMPERATURE), 0.0)
        self.coarse.ProcessInfo[Kratos.STEP] = 2
        process.ExecuteFinalizeSolutionStep()
        self.assertNotEqual(self.fine.GetNode(1).GetSolutionStepValue(Kratos.TEMPERATURE), 0.0)

    def test_InvalidBoundingBoxRejected(self):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "coarse_model_part_name" : "Coarse",
                "fine_model_part_name"   : "Fine",
                "input_fields"           : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"          : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "bounding_box"           : [0.0, 0.0, 0.0]
            }
        }""")
        with self.assertRaisesRegex(ValueError, "bounding_box"):
            superresolution_process.Factory(settings, self.model)


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestGridOperators2D(KratosUnittest.TestCase):
    """The thin-axis (squeeze) idiom and the modafno interface on
    GridInferenceProcess/SuperResolutionProcess."""

    def setUp(self):
        self.checkpoint = Path("test_grid_operators_2d.pt")
        self.model = Kratos.Model()
        self.part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=3,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        for node in self.part.Nodes:  # z-independent linear field
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0 + 2.0 * node.X + 3.0 * node.Y)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def _CreateProcess(self, model_interface="grid"):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import grid_inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : { "checkpoint_file" : "%s", "device" : "cpu" },
                "model_interface" : "%s",
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "grid_shape"      : [6, 6, 2],
                "squeeze_axis"    : 2
            }
        }""" % (self.checkpoint, model_interface))
        return grid_inference_process.Factory(settings, self.model)

    def test_SqueezeAxisIdentityExact(self):
        class Identity2D(torch.nn.Module):
            def forward(self, x):  # (1, C, A, B) -> (1, C, A, B)
                return x

        torch.jit.script(Identity2D()).save(str(self.checkpoint))
        process = self._CreateProcess()
        self.part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        for node in self.part.Nodes:
            expected = 1.0 + 2.0 * node.X + 3.0 * node.Y
            self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.TEMPERATURE), expected, places=8)

    def test_ModAfnoInterfacePassesTime(self):
        class AddMod(torch.nn.Module):
            def forward(self, x, mod):  # ModAFNO's (x, mod) contract
                return x + mod[0, 0]

        torch.jit.script(AddMod()).save(str(self.checkpoint))
        process = self._CreateProcess(model_interface="modafno")
        self.part.ProcessInfo[Kratos.STEP] = 1
        self.part.ProcessInfo[Kratos.TIME] = 2.5
        process.ExecuteFinalizeSolutionStep()

        for node in self.part.Nodes:
            expected = 1.0 + 2.0 * node.X + 3.0 * node.Y + 2.5
            self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.TEMPERATURE), expected, places=6)

    def test_UnknownInterfaceRaises(self):
        class Identity2D(torch.nn.Module):
            def forward(self, x):
                return x

        torch.jit.script(Identity2D()).save(str(self.checkpoint))
        with self.assertRaisesRegex(ValueError, "model interface"):
            self._CreateProcess(model_interface="voxel")


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestGridOperatorZoo(KratosUnittest.TestCase):
    """Real 2D grid operators (FNO dimension=2, AFNO) through the squeeze
    path of GridInferenceProcess."""

    def setUp(self):
        self.checkpoint = Path("test_grid_operator_zoo.mdlus")
        self.model = Kratos.Model()
        self.part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=3,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        for node in self.part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + node.Y)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def _RunThroughProcess(self, grid_shape="[8, 8, 2]"):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import grid_inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "grid_shape"      : %s,
                "squeeze_axis"    : 2
            }
        }""" % (self.checkpoint, grid_shape))
        process = grid_inference_process.Factory(settings, self.model)
        self.part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        values = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.part.Nodes])
        self.assertTrue(numpy.isfinite(values).all())
        self.assertGreater(numpy.abs(values).max(), 0.0)

    def test_Fno2DThroughProcess(self):
        from physicsnemo.models.fno import FNO

        torch.manual_seed(0)
        fno = FNO(in_channels=1, out_channels=1, dimension=2,
                  latent_channels=8, num_fno_layers=2, num_fno_modes=2, padding=2)
        fno.save(str(self.checkpoint))
        self._RunThroughProcess()

    def test_DlwpThroughProcess(self):
        """DLWP's five-dimensional (B, C, 6, H, H) cubed-sphere layout maps
        onto the 3D (C, D, H, W) path with D as the face axis - here on a
        [6, 8, 8] grid without the squeeze idiom (physics semantics aside:
        the cubed-sphere padding stays weather-specific)."""
        from physicsnemo.models.dlwp import DLWP
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import grid_inference_process

        torch.manual_seed(0)
        DLWP(nr_input_channels=1, nr_output_channels=1, nr_initial_channels=4, depth=1).save(
            str(self.checkpoint))
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "grid_shape"      : [6, 8, 8]
            }
        }""" % self.checkpoint)
        process = grid_inference_process.Factory(settings, self.model)
        self.part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        values = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.part.Nodes])
        self.assertTrue(numpy.isfinite(values).all())
        self.assertGreater(numpy.abs(values).max(), 0.0)

    def test_AfnoThroughProcess(self):
        from physicsnemo.models.afno import AFNO

        torch.manual_seed(0)
        afno = AFNO(inp_shape=[8, 8], in_channels=1, out_channels=1,
                    patch_size=[4, 4], embed_dim=16, depth=2, num_blocks=2)
        afno.save(str(self.checkpoint))
        self._RunThroughProcess()


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestSRResNetThroughProcess(KratosUnittest.TestCase):
    """Runs a real (untrained) SRResNet through the process and exercises the
    physicsnemo checkpoint path of model_registry."""

    def setUp(self):
        self.checkpoint = Path("test_srresnet.mdlus")
        model = SRResNet(in_channels=1, out_channels=1, conv_layer_size=4,
                         n_resid_blocks=1, scaling_factor=2)
        model.save(str(self.checkpoint))

        self.model = Kratos.Model()
        self.coarse = CreateStructuredTetModelPart(self.model, "Coarse", divisions=3)
        self.fine = CreateStructuredTetModelPart(self.model, "Fine", divisions=5,
                                                 historical_variables=(Kratos.TEMPERATURE,))
        for node in self.coarse.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def test_PhysicsNemoCheckpointLoads(self):
        settings = Kratos.Parameters("""{
            "checkpoint_file" : "test_srresnet.mdlus",
            "checkpoint_type" : "physicsnemo",
            "device"          : "cpu"
        }""")
        model, device = model_registry.LoadModel(settings)
        out = model(torch.zeros(1, 1, 4, 4, 4))
        self.assertEqual(tuple(out.shape), (1, 1, 8, 8, 8))

    def test_SRResNetThroughProcess(self):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "coarse_model_part_name" : "Coarse",
                "fine_model_part_name"   : "Fine",
                "model_settings"         : {
                    "checkpoint_file" : "test_srresnet.mdlus",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "input_fields"           : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"          : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "coarse_grid_shape"      : [6, 6, 6]
            }
        }""")
        process = superresolution_process.Factory(settings, self.model)
        self.coarse.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()  # untrained net: only plumbing/shapes are checked
        values = [node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.fine.Nodes]
        self.assertTrue(numpy.isfinite(values).all())


if __name__ == '__main__':
    KratosUnittest.main()
