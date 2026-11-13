"""Tests for two architecture recipes: DPOT through SequenceInferenceProcess
and TopoDiff through the diffusion bridge."""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.models.dpot import DPOTNet
    have_dpot = True
except ImportError:
    have_dpot = False

try:
    from physicsnemo.models.topodiff import TopoDiff
    have_topodiff = True
except ImportError:
    have_topodiff = False


def _TinyDpot(in_timesteps=2, out_timesteps=2, seed=0):
    torch.manual_seed(seed)
    return DPOTNet(inp_shape=8, patch_size=2, in_channels=1, out_channels=1,
                   in_timesteps=in_timesteps, out_timesteps=out_timesteps,
                   embed_dim=16, depth=1, num_blocks=1, modes=4)


def _TinyTopoDiff(seed=0):
    torch.manual_seed(seed)
    # model_channels below 64 makes the attention layer derive zero heads
    return TopoDiff(img_resolution=16, in_channels=4, out_channels=1, model_channels=64,
                    channel_mult=[1, 1], num_blocks=1, attn_resolutions=[])


@KratosUnittest.skipUnless(have_torch and have_dpot,
                           "Missing required python modules: torch, physicsnemo.")
class TestDpotInterface(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_dpot.mdlus")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=3,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.X + node.Y)

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting(str(self.checkpoint))

    def _Process(self, window_size=2, checkpoint=None):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import (
            sequence_inference_process)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"     : "Main",
                "model_settings"      : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "model_interface"     : "dpot",
                "input_fields"        : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"       : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "grid_shape"          : [8, 8, 2],
                "squeeze_axis"        : 2,
                "window_as_time_axis" : true,
                "window_size"         : %d
            }
        }""" % (checkpoint or self.checkpoint, window_size))
        return sequence_inference_process.Factory(settings, self.model)

    def test_TheRolloutWritesFiniteStates(self):
        _TinyDpot().save(str(self.checkpoint))
        process = self._Process()
        for step in (1, 2, 3):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
        written = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        self.assertTrue(numpy.isfinite(written).all())
        self.assertGreater(numpy.abs(written).max(), 0.0)

    def test_TheAxisPermutationIsTheWholeAdapter(self):
        """DPOT reads (B, H, W, T, C) where every other grid model here reads
        (B, C, T, H, W). A model that returns its input unchanged must come
        back as the sampled field itself - so the two permutations undo each
        other rather than quietly transposing the domain."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge

        class Identity(torch.nn.Module):
            in_timesteps = 2

            def forward(self, x):       # (B, H, W, T, C) in and out
                return x

        process = self._Process()
        process._model = Identity()
        process._device = torch.device("cpu")
        process._normalization = None
        for step in (1, 2, 3):          # two to fill the window, one to write
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()

        # the reference is the field's own GRID ROUND TRIP: sampling onto the
        # lattice and scattering back is interpolation, not the identity on
        # nodal values, so comparing against the raw nodal field would be wrong
        bounding_box = grid_bridge.ComputeBoundingBox(self.model_part)
        sampled, _ = grid_bridge.SampleFieldsOnGrid(
            self.model_part, [("PRESSURE", "node_historical")], (8, 8, 2), bounding_box)
        squeezed = sampled.mean(axis=3)
        expected_grid = numpy.repeat(numpy.expand_dims(squeezed, 3), 2, axis=3)
        for node in self.model_part.Nodes:
            node.SetValue(Kratos.NODAL_PAUX, 0.0)
        grid_bridge.ScatterGridToNodes(
            expected_grid, bounding_box, self.model_part,
            [("NODAL_PAUX", "node_non_historical")])

        written = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        expected = numpy.array([
            node.GetValue(Kratos.NODAL_PAUX) for node in self.model_part.Nodes])
        numpy.testing.assert_allclose(written, expected, atol=1e-9)

    def test_AWindowThatDisagreesWithTheModelIsReported(self):
        _TinyDpot(in_timesteps=3).save(str(self.checkpoint))
        process = self._Process(window_size=2)
        with self.assertRaisesRegex(ValueError, "in_timesteps"):
            for step in (1, 2):
                self.model_part.ProcessInfo[Kratos.STEP] = step
                process.ExecuteFinalizeSolutionStep()

    def test_AVolumetricGridIsRefused(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import (
            sequence_inference_process)

        _TinyDpot().save(str(self.checkpoint))
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"     : "Main",
                "model_settings"      : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "model_interface"     : "dpot",
                "input_fields"        : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"       : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "grid_shape"          : [8, 8, 8],
                "window_as_time_axis" : true,
                "window_size"         : 2
            }
        }""" % self.checkpoint)
        process = sequence_inference_process.Factory(settings, self.model)
        with self.assertRaisesRegex(ValueError, "2-D grid series"):
            for step in (1, 2):
                self.model_part.ProcessInfo[Kratos.STEP] = step
                process.ExecuteFinalizeSolutionStep()

    def test_AnUnknownInterfaceRaises(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import (
            sequence_inference_process)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_interface" : "fourier",
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""")
        with self.assertRaisesRegex(ValueError, "model interface"):
            sequence_inference_process.Factory(settings, self.model)


@KratosUnittest.skipUnless(have_torch and have_topodiff,
                           "Missing required python modules: torch, physicsnemo.")
class TestTopoDiffDenoiser(KratosUnittest.TestCase):
    def test_TheAdapterPassesConstraintsNatively(self):
        """TopoDiff takes its constraint channels as a second ARGUMENT,
        where DiT wants them concatenated into the input."""
        from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils

        model = _TinyTopoDiff().eval()
        wrapped = diffusion_utils.WrapDiffusionModel(model, "topodiff", out_channels=1)
        x = torch.randn(2, 1, 16, 16)
        constraints = torch.randn(2, 3, 16, 16)
        t = torch.full((2,), 0.4)
        with torch.no_grad():
            reference = model(x, constraints, t)
            through_adapter = wrapped(x, t, condition=constraints)
        torch.testing.assert_close(through_adapter, reference)

    def test_TheChannelCountMustBeExplicit(self):
        """TopoDiff exposes no out_channels attribute, so the bridge cannot
        infer it and says so rather than guessing."""
        from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils

        model = _TinyTopoDiff()
        self.assertFalse(hasattr(model, "out_channels"))
        with self.assertRaisesRegex(ValueError, "out_channels"):
            diffusion_utils.WrapDiffusionModel(model, "topodiff")

    def test_AMissingConditionIsReported(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils

        wrapped = diffusion_utils.WrapDiffusionModel(
            _TinyTopoDiff().eval(), "topodiff", out_channels=1)
        with self.assertRaisesRegex(ValueError, "condition"):
            wrapped(torch.randn(1, 1, 16, 16), torch.full((1,), 0.5))

    def test_ItGeneratesAnEnsembleThroughTheSampler(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils

        wrapped = diffusion_utils.WrapDiffusionModel(
            _TinyTopoDiff().eval(), "topodiff", out_channels=1)
        ensemble = diffusion_utils.GenerateEnsemble(
            wrapped, numpy.zeros((3, 16, 16)), Kratos.Parameters("""{
                "num_samples"        : 2,
                "num_steps"          : 4,
                "seed"               : 0,
                "denoiser_interface" : "protocol",
                "output_channels"    : 1
            }"""))
        self.assertEqual(ensemble.shape, (2, 1, 16, 16))
        self.assertTrue(numpy.isfinite(ensemble).all())

    def test_TooFewModelChannelsIsAnUpstreamTrap(self):
        """Worth pinning because the message names num_heads, not width:
        TopoDiff derives its head count from model_channels, so anything
        below 64 builds zero heads and fails at construction."""
        torch.manual_seed(0)
        with self.assertRaisesRegex(ValueError, "num_heads"):
            TopoDiff(img_resolution=16, in_channels=4, out_channels=1, model_channels=32,
                     channel_mult=[1, 1], num_blocks=1, attn_resolutions=[])


if __name__ == '__main__':
    KratosUnittest.main()
