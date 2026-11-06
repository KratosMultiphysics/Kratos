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
    import physicsnemo.models.rnn.rnn_one2many
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


def _WriteGridSeries(directory, steps, shape=(1, 2, 2, 2)):
    directory.mkdir(parents=True, exist_ok=True)
    for step in steps:
        numpy.savez(
            directory / f"grid_{step}.npz",
            grid=numpy.full(shape, float(step), dtype=numpy.float32),
            TIME=float(step), STEP=step,
            bounding_box=numpy.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0]))


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestGridSequenceDataset(KratosUnittest.TestCase):
    def setUp(self):
        self.directory = Path("test_grid_sequence_dataset")
        _WriteGridSeries(self.directory, range(5))

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.directory))

    def test_SequenceLayout(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import CreateGridSequenceDataset

        dataset = CreateGridSequenceDataset(self.directory, nr_tsteps=2)
        self.assertEqual(len(dataset), 3)
        initial, future = dataset[1]
        self.assertEqual(tuple(initial.shape), (1, 1, 2, 2, 2))
        self.assertEqual(tuple(future.shape), (1, 2, 2, 2, 2))
        self.assertTrue(bool((initial == 1.0).all()))
        self.assertTrue(bool((future[:, 0] == 2.0).all()))
        self.assertTrue(bool((future[:, 1] == 3.0).all()))

    def test_SqueezeAxis(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import CreateGridSequenceDataset

        dataset = CreateGridSequenceDataset(self.directory, nr_tsteps=1, squeeze_axis=2)
        initial, future = dataset[0]
        self.assertEqual(tuple(initial.shape), (1, 1, 2, 2))
        self.assertEqual(tuple(future.shape), (1, 1, 2, 2))

    def test_TooFewGridsRaise(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import CreateGridSequenceDataset

        with self.assertRaisesRegex(ValueError, "nr_tsteps"):
            CreateGridSequenceDataset(self.directory, nr_tsteps=5)
        with self.assertRaisesRegex(ValueError, "squeeze_axis"):
            CreateGridSequenceDataset(self.directory, nr_tsteps=1, squeeze_axis=3)


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestSequenceInferenceProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_sequence_toy_model.pt")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 2.0)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def _SaveToyModel(self):
        class Roll(torch.nn.Module):
            def forward(self, x):  # (N, C, 1, D, H, W) -> (N, C, 3, D, H, W)
                return torch.cat([x * 1.0, x * 2.0, x * 3.0], dim=2)

        torch.jit.script(Roll()).save(str(self.checkpoint))

    def _CreateProcess(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import sequence_inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : { "checkpoint_file" : "%s", "device" : "cpu" },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "grid_shape"      : [4, 4, 4]
            }
        }""" % self.checkpoint)
        return sequence_inference_process.Factory(settings, self.model)

    def test_SeedThenRollout(self):
        self._SaveToyModel()
        process = self._CreateProcess()

        def temperature(node_id=1):
            return self.model_part.GetNode(node_id).GetSolutionStepValue(Kratos.TEMPERATURE)

        # step 1 seeds the rollout: 3 states buffered, nothing written yet
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        self.assertEqual(process.predicted_steps_left, 3)
        self.assertEqual(temperature(), 0.0)

        # steps 2..4 pop the buffered states: PRESSURE * (1, 2, 3)
        for step, expected in ((2, 2.0), (3, 4.0), (4, 6.0)):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
            self.assertAlmostEqual(temperature(), expected, places=10)

        # step 5: exhausted - warns once, leaves the last state in place
        self.model_part.ProcessInfo[Kratos.STEP] = 5
        process.ExecuteFinalizeSolutionStep()
        self.assertAlmostEqual(temperature(), 6.0, places=10)
        self.assertEqual(process.predicted_steps_left, 0)

    def test_OutputNormalizationFromTheCardIsApplied(self):
        """Every buffered state is de-normalized on its way out - the card
        belongs to the rollout, not just to the seeding step."""
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
        self._SaveToyModel()
        mean, std = 1.5, 10.0
        model_registry.SaveModelCard(self.checkpoint, {
            "output_normalization": {"type": "mean_std", "mean": [mean], "std": [std]}})
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")
        process = self._CreateProcess()

        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()  # seeds: PRESSURE * (1, 2, 3) buffered
        for step, raw in ((2, 2.0), (3, 4.0), (4, 6.0)):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
            self.assertAlmostEqual(
                self.model_part.GetNode(1).GetSolutionStepValue(Kratos.TEMPERATURE),
                std * raw + mean, places=10)

    def test_InputNormalizationFromTheCardIsApplied(self):
        # sampled grids are standardized before they seed the rollout:
        # states = ((P - mean) / std) * (1, 2, 3)
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
        self._SaveToyModel()
        mean, std = 1.0, 4.0
        model_registry.SaveModelCard(self.checkpoint, {
            "input_normalization": {"type": "mean_std", "mean": [mean], "std": [std]}})
        self.addCleanup(KratosUtilities.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")
        process = self._CreateProcess()

        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        for step, factor in ((2, 1.0), (3, 2.0), (4, 3.0)):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
            self.assertAlmostEqual(
                self.model_part.GetNode(1).GetSolutionStepValue(Kratos.TEMPERATURE),
                (2.0 - mean) / std * factor, places=10)

    def test_WrongOutputRankRaises(self):
        class Flat(torch.nn.Module):
            def forward(self, x):
                return x[:, :, 0]

        torch.jit.script(Flat()).save(str(self.checkpoint))
        process = self._CreateProcess()
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        with self.assertRaisesRegex(ValueError, "sequence per sample"):
            process.ExecuteFinalizeSolutionStep()


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestWindowAsTimeAxis(KratosUnittest.TestCase):
    """window_as_time_axis: the sampled window becomes the model's time axis
    (spatiotemporal block operators, FNO dimension=4)."""

    def setUp(self):
        self.checkpoint = Path("test_sequence_window_model.pt")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def _SetPressure(self, value):
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, value)

    def _CreateProcess(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import sequence_inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"     : "Main",
                "model_settings"      : { "checkpoint_file" : "%s", "device" : "cpu" },
                "input_fields"        : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"       : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "grid_shape"          : [4, 4, 4],
                "window_as_time_axis" : true,
                "window_size"         : 2
            }
        }""" % self.checkpoint)
        return sequence_inference_process.Factory(settings, self.model)

    def test_WindowSeedAndRollout(self):
        class Identity(torch.nn.Module):
            def forward(self, x):  # (1, C, K, D, H, W) -> same block back
                return x

        torch.jit.script(Identity()).save(str(self.checkpoint))
        process = self._CreateProcess()

        def temperature():
            return self.model_part.GetNode(1).GetSolutionStepValue(Kratos.TEMPERATURE)

        # step 1: first window entry, no prediction
        self._SetPressure(10.0)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        self.assertEqual(process.predicted_steps_left, 0)
        self.assertEqual(temperature(), 0.0)

        # step 2: window full -> seeded with the (10, 20) block
        self._SetPressure(20.0)
        self.model_part.ProcessInfo[Kratos.STEP] = 2
        process.ExecuteFinalizeSolutionStep()
        self.assertEqual(process.predicted_steps_left, 2)

        # steps 3, 4 pop the identity-predicted block: 10 then 20
        for step, expected in ((3, 10.0), (4, 20.0)):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
            self.assertAlmostEqual(temperature(), expected, places=10)

    def test_WindowSizeValidation(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import sequence_inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"     : "Main",
                "input_fields"        : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"       : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                "window_as_time_axis" : true,
                "window_size"         : 1
            }
        }""")
        with self.assertRaisesRegex(ValueError, "window_size"):
            sequence_inference_process.Factory(settings, self.model)

    @KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
    def test_Fno4DThroughProcess(self):
        from physicsnemo.models.fno import FNO

        torch.manual_seed(0)
        fno = FNO(in_channels=1, out_channels=1, dimension=4,
                  latent_channels=4, num_fno_layers=2, num_fno_modes=[1, 2, 2, 2], padding=0)
        checkpoint = Path("test_sequence_fno4d.mdlus")
        fno.save(str(checkpoint))
        try:
            from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import sequence_inference_process
            settings = Kratos.Parameters("""{
                "Parameters": {
                    "model_part_name"     : "Main",
                    "model_settings"      : {
                        "checkpoint_file" : "test_sequence_fno4d.mdlus",
                        "checkpoint_type" : "physicsnemo",
                        "device"          : "cpu"
                    },
                    "input_fields"        : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                    "output_fields"       : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                    "grid_shape"          : [4, 4, 4],
                    "window_as_time_axis" : true,
                    "window_size"         : 2
                }
            }""")
            process = sequence_inference_process.Factory(settings, self.model)

            self._SetPressure(1.0)
            self.model_part.ProcessInfo[Kratos.STEP] = 1
            process.ExecuteFinalizeSolutionStep()
            self._SetPressure(2.0)
            self.model_part.ProcessInfo[Kratos.STEP] = 2
            process.ExecuteFinalizeSolutionStep()
            # FNO(dimension=4) preserves the block: K predicted states buffered
            self.assertEqual(process.predicted_steps_left, 2)

            self.model_part.ProcessInfo[Kratos.STEP] = 3
            process.ExecuteFinalizeSolutionStep()
            values = numpy.array([
                node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
            self.assertTrue(numpy.isfinite(values).all())
            self.assertGreater(numpy.abs(values).max(), 0.0)
        finally:
            KratosUtilities.DeleteFileIfExisting("test_sequence_fno4d.mdlus")


@KratosUnittest.skipUnless(have_torch and have_physicsnemo,
                           "Missing required python modules: torch, physicsnemo.")
class TestOne2ManyRNNThroughProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_sequence_rnn_model.mdlus")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0 + node.X)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def test_RealRNNRollout(self):
        from physicsnemo.models.rnn.rnn_one2many import One2ManyRNN
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import sequence_inference_process
        torch.manual_seed(0)
        rnn = One2ManyRNN(
            input_channels=1, dimension=3, nr_latent_channels=4,
            nr_residual_blocks=1, nr_downsamples=1, nr_tsteps=2)
        rnn.save(str(self.checkpoint))

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
                "grid_shape"      : [8, 8, 8]
            }
        }""" % self.checkpoint)
        process = sequence_inference_process.Factory(settings, self.model)

        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()  # seed
        self.assertEqual(process.predicted_steps_left, 2)

        self.model_part.ProcessInfo[Kratos.STEP] = 2
        process.ExecuteFinalizeSolutionStep()  # first predicted state
        values = numpy.array([
            node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in self.model_part.Nodes])
        self.assertTrue(numpy.isfinite(values).all())
        self.assertGreater(numpy.abs(values).max(), 0.0)


if __name__ == '__main__':
    KratosUnittest.main()
