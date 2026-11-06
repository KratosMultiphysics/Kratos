"""Tests for training streamed out of a running solve, and for shrink-and-
perturb warm restarts."""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.processes.export import dataset_export_process
from KratosMultiphysics.PhysicsNeMoApplication.training import streaming_dataset
from KratosMultiphysics.PhysicsNeMoApplication.processes.export import (
    streaming_dataset_export_process)
from KratosMultiphysics.PhysicsNeMoApplication.training import torch_dataset
from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
sys.path.insert(0, str(Path(__file__).parent / "kratos_solver_cases"))

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.datapipes import IterableDatasetBase  # noqa: F401
    have_iterable = have_torch
except ImportError:
    have_iterable = False

try:
    from physicsnemo.nn import shrink_and_perturb_  # noqa: F401
    have_shrink_perturb = have_torch
except ImportError:
    have_shrink_perturb = False

_MISSING = "Missing required python modules: torch, physicsnemo >= 2.2."


def _CreateModelPart(model, name="Stream", n=6):
    model_part = model.CreateModelPart(name)
    model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
    model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
    model_part.SetBufferSize(1)
    for i in range(n):
        model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
    return model_part


def _StepValues(model_part, step):
    for node in model_part.Nodes:
        node.SetSolutionStepValue(Kratos.PRESSURE, float(step) + node.X)
        node.SetSolutionStepValue(Kratos.TEMPERATURE, 2.0 * (float(step) + node.X))


class TestLiveSampleQueue(KratosUnittest.TestCase):
    """Pure python: no torch or physicsnemo needed."""

    def test_DrainYieldsInOrderAndStopsAtClose(self):
        queue = streaming_dataset.LiveSampleQueue()
        for i in range(3):
            queue.Push({"STEP": numpy.array(i)})
        queue.Close()
        drained = [int(sample["STEP"]) for sample in queue.Drain()]
        self.assertEqual(drained, [0, 1, 2])

    def test_CloseIsRequiredForTerminationAndIsOneWay(self):
        queue = streaming_dataset.LiveSampleQueue()
        self.assertFalse(queue.closed)
        queue.Push({"STEP": numpy.array(0)})
        queue.Close()
        self.assertTrue(queue.closed)
        # a closed queue refuses more work rather than silently dropping it
        with self.assertRaisesRegex(RuntimeError, "Push after Close"):
            queue.Push({"STEP": numpy.array(1)})

    def test_Validation(self):
        with self.assertRaisesRegex(ValueError, "max_size"):
            streaming_dataset.LiveSampleQueue(-1)


@KratosUnittest.skipUnless(have_iterable, _MISSING)
class TestStreamingDataset(KratosUnittest.TestCase):

    def _Fill(self, queue, steps=3):
        model = Kratos.Model()
        model_part = _CreateModelPart(model)
        for step in range(1, steps + 1):
            model_part.ProcessInfo[Kratos.STEP] = step
            model_part.ProcessInfo[Kratos.TIME] = 0.1 * step
            _StepValues(model_part, step)
            queue.Push(streaming_dataset.GatherSampleArrays(
                model_part, [("PRESSURE", "node_historical"),
                             ("TEMPERATURE", "node_historical")]))
        queue.Close()

    def test_ItemsMatchTheNpzDatasetLayout(self):
        queue = streaming_dataset.LiveSampleQueue()
        self._Fill(queue)
        dataset = streaming_dataset.CreateStreamingDataset(
            queue, ["PRESSURE__node_historical"], ["TEMPERATURE__node_historical"])

        items = list(iter(dataset))
        self.assertEqual(len(items), 3)
        inputs, targets = items[0]
        self.assertEqual(tuple(inputs.shape), (6, 1))
        self.assertEqual(tuple(targets.shape), (6, 1))
        self.assertEqual(inputs.dtype, torch.float32)
        # targets are exactly 2x the inputs by construction
        self.assertLess(float((targets - 2.0 * inputs).abs().max()), 1e-6)

    def test_IsBothAPhysicsNemoAndATorchIterableDataset(self):
        """The dual inheritance is load-bearing: a bare IterableDatasetBase
        subclass is rejected by torch's DataLoader for having no len()."""
        queue = streaming_dataset.LiveSampleQueue()
        self._Fill(queue, steps=2)
        dataset = streaming_dataset.CreateStreamingDataset(
            queue, ["PRESSURE__node_historical"], ["TEMPERATURE__node_historical"])
        self.assertIsInstance(dataset, IterableDatasetBase)
        self.assertIsInstance(dataset, torch.utils.data.IterableDataset)
        self.assertTrue(dataset.yields_batches)

        loader = torch.utils.data.DataLoader(dataset, batch_size=None)
        batches = list(loader)
        self.assertEqual(len(batches), 2)
        self.assertEqual(tuple(batches[0][0].shape), (6, 1))

    def test_TargetsSurvivePhysicsNemoCollation(self):
        """yields_batches=True bypasses the (data, metadata) unpacking that
        would otherwise silently discard the targets."""
        from physicsnemo.datapipes import DataLoader as PhysicsNemoDataLoader

        queue = streaming_dataset.LiveSampleQueue()
        self._Fill(queue, steps=2)
        dataset = streaming_dataset.CreateStreamingDataset(
            queue, ["PRESSURE__node_historical"], ["TEMPERATURE__node_historical"])
        items = list(PhysicsNemoDataLoader(dataset, batch_size=4))
        self.assertEqual(len(items), 2)
        inputs, targets = items[0]
        self.assertEqual(tuple(inputs.shape), (6, 1))
        self.assertEqual(tuple(targets.shape), (6, 1))

    def test_WorkerDuplicationRefused(self):
        queue = streaming_dataset.LiveSampleQueue()
        self._Fill(queue, steps=2)
        dataset = streaming_dataset.CreateStreamingDataset(
            queue, ["PRESSURE__node_historical"], ["TEMPERATURE__node_historical"])
        loader = torch.utils.data.DataLoader(dataset, batch_size=None, num_workers=2)
        with self.assertRaises(RuntimeError):
            list(loader)

    def test_Validation(self):
        queue = streaming_dataset.LiveSampleQueue()
        with self.assertRaisesRegex(ValueError, "non-empty"):
            streaming_dataset.CreateStreamingDataset(queue, [], ["A__node_historical"])
        queue.Push({"STEP": numpy.array(1)})
        queue.Close()
        dataset = streaming_dataset.CreateStreamingDataset(
            queue, ["MISSING__node_historical"], ["MISSING__node_historical"])
        with self.assertRaisesRegex(KeyError, "MISSING"):
            list(iter(dataset))


@KratosUnittest.skipUnless(have_iterable, _MISSING)
class TestStreamingMatchesTheFileRoundTrip(KratosUnittest.TestCase):
    """The point of the whole feature: streaming and dumping must agree."""

    def setUp(self):
        self.output_path = Path("test_streaming_dump")

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))

    def test_SameSolveYieldsIdenticalSamplesEitherWay(self):
        fields = """[ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" },
                      { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]"""

        # (a) the shipped path: dump .npz per step, read them back
        model = Kratos.Model()
        model_part = _CreateModelPart(model, "Dumped")
        dumper = dataset_export_process.Factory(Kratos.Parameters("""{
            "Parameters": { "model_part_name" : "Dumped",
                            "list_of_fields"  : %s,
                            "output_path"     : "%s" }
        }""" % (fields, self.output_path)), model)
        dumper.ExecuteInitialize()

        # (b) the streaming path: push into a queue
        streamed_model = Kratos.Model()
        streamed_part = _CreateModelPart(streamed_model, "Streamed")
        streamer = streaming_dataset_export_process.Factory(Kratos.Parameters("""{
            "Parameters": { "model_part_name" : "Streamed",
                            "list_of_fields"  : %s }
        }""" % fields), streamed_model)

        for step in range(1, 5):
            for part in (model_part, streamed_part):
                part.ProcessInfo[Kratos.STEP] = step
                part.ProcessInfo[Kratos.TIME] = 0.1 * step
                _StepValues(part, step)
            dumper.ExecuteFinalizeSolutionStep()
            streamer.ExecuteFinalizeSolutionStep()
        streamer.ExecuteFinalize()

        keys_in = ["PRESSURE__node_historical"]
        keys_out = ["TEMPERATURE__node_historical"]
        dumped = torch_dataset.CreateNpzDataset(self.output_path, keys_in, keys_out)
        streamed = list(iter(streaming_dataset.CreateStreamingDataset(
            streamer.queue, keys_in, keys_out)))

        self.assertEqual(len(streamed), len(dumped))
        for index, (inputs, targets) in enumerate(streamed):
            reference_inputs, reference_targets = dumped[index]
            self.assertTrue(torch.equal(inputs, reference_inputs),
                            msg=f"inputs differ at sample {index}")
            self.assertTrue(torch.equal(targets, reference_targets),
                            msg=f"targets differ at sample {index}")

    def test_TrainingConsumesTheStream(self):
        model = Kratos.Model()
        model_part = _CreateModelPart(model, "TrainStream")
        streamer = streaming_dataset_export_process.Factory(Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "TrainStream",
                "list_of_fields"  : [
                    { "variable_name" : "PRESSURE",    "data_location" : "node_historical" },
                    { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }"""), model)
        for step in range(1, 9):
            model_part.ProcessInfo[Kratos.STEP] = step
            model_part.ProcessInfo[Kratos.TIME] = 0.1 * step
            _StepValues(model_part, step)
            streamer.ExecuteFinalizeSolutionStep()
        streamer.ExecuteFinalize()

        dataset = streaming_dataset.CreateStreamingDataset(
            streamer.queue, ["PRESSURE__node_historical"], ["TEMPERATURE__node_historical"])
        torch.manual_seed(0)
        surrogate = torch.nn.Linear(1, 1)
        history = training_utils.TrainModel(surrogate, dataset, Kratos.Parameters("""{
            "epochs"    : 1,
            "streaming" : true,
            "shuffle"   : false,
            "device"    : "cpu",
            "seed"      : 0
        }"""))
        self.assertEqual(len(history), 1)
        self.assertTrue(numpy.isfinite(history[0]))
        self.assertEqual(dataset.emitted, 8)

    def test_StreamingRejectsShuffleAndExtraEpochs(self):
        queue = streaming_dataset.LiveSampleQueue()
        queue.Close()
        dataset = streaming_dataset.CreateStreamingDataset(
            queue, ["PRESSURE__node_historical"], ["TEMPERATURE__node_historical"])
        model = torch.nn.Linear(1, 1)
        with self.assertRaisesRegex(ValueError, "shuffle"):
            training_utils.TrainModel(model, dataset, Kratos.Parameters(
                '{"epochs": 1, "streaming": true, "shuffle": true, "device": "cpu"}'))
        with self.assertRaisesRegex(ValueError, "epochs"):
            training_utils.TrainModel(model, dataset, Kratos.Parameters(
                '{"epochs": 3, "streaming": true, "shuffle": false, "device": "cpu"}'))


@KratosUnittest.skipUnless(have_shrink_perturb, _MISSING)
class TestWarmRestart(KratosUnittest.TestCase):

    def setUp(self):
        self.checkpoint = "test_warm_restart.pt"
        torch.manual_seed(0)
        self.dataset = torch.utils.data.TensorDataset(
            torch.randn(64, 3), torch.randn(64, 1))

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(self.checkpoint)

    def _Model(self):
        torch.manual_seed(1)
        return torch.nn.Sequential(torch.nn.Linear(3, 8), torch.nn.LayerNorm(8),
                                   torch.nn.Tanh(), torch.nn.Linear(8, 1))

    def test_IdentityAtShrinkOneAndZeroPerturb(self):
        model = self._Model()
        before = [p.clone() for p in model.parameters()]
        training_utils.TrainModel(model, self.dataset, Kratos.Parameters("""{
            "epochs" : 0, "device" : "cpu",
            "warm_restart" : { "shrink" : 1.0, "perturb" : 0.0 }
        }"""))
        for original, current in zip(before, model.parameters()):
            self.assertTrue(torch.equal(original, current))

    def test_ShrinkMovesWeightsTowardInitialization(self):
        model = self._Model()
        weight_before = model[0].weight.clone()
        training_utils.TrainModel(model, self.dataset, Kratos.Parameters("""{
            "epochs" : 0, "device" : "cpu", "seed" : 7,
            "warm_restart" : { "shrink" : 0.5, "perturb" : 0.0 }
        }"""))
        self.assertLess(float(model[0].weight.norm()), float(weight_before.norm()))
        self.assertAlmostEqual(
            float((model[0].weight - 0.5 * weight_before).abs().max()), 0.0, places=12)

    def test_ReproducibleUnderAFixedSeed(self):
        results = []
        for _ in range(2):
            model = self._Model()
            training_utils.TrainModel(model, self.dataset, Kratos.Parameters("""{
                "epochs" : 0, "device" : "cpu", "seed" : 11,
                "warm_restart" : { "shrink" : 0.5, "perturb" : 0.1 }
            }"""))
            results.append(model[0].weight.clone())
        self.assertTrue(torch.equal(results[0], results[1]))

    def test_NormLayersAndBiasesAreLeftAloneByDefault(self):
        """Upstream perturbs everything, which halves LayerNorm gains toward
        zero and crashes outright on integer parameters."""
        model = self._Model()
        norm_weight = model[1].weight.clone()
        bias = model[0].bias.clone()
        training_utils.TrainModel(model, self.dataset, Kratos.Parameters("""{
            "epochs" : 0, "device" : "cpu", "seed" : 3,
            "warm_restart" : { "shrink" : 0.5, "perturb" : 0.1 }
        }"""))
        self.assertTrue(torch.equal(model[1].weight, norm_weight))
        self.assertTrue(torch.equal(model[0].bias, bias))
        # ... while the 2-D weights did move
        self.assertFalse(torch.equal(model[0].weight, self._Model()[0].weight))

    def test_Validation(self):
        model = self._Model()
        for block, pattern in (
                ('{"shrink": 1.5}', "shrink"),
                ('{"shrink": -0.1}', "shrink"),
                ('{"noise": "uniform", "shrink": 0.5}', "noise")):
            with self.assertRaisesRegex(ValueError, pattern):
                training_utils.TrainModel(model, self.dataset, Kratos.Parameters(
                    '{"epochs": 0, "device": "cpu", "warm_restart": %s}' % block))

    def test_RoundTripThroughACheckpoint(self):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
        model = self._Model()
        first = training_utils.TrainModel(model, self.dataset, Kratos.Parameters(
            '{"epochs": 5, "device": "cpu", "seed": 0}'))
        training_utils.SaveTrainedModel(model, self.checkpoint)

        restored, _ = model_registry.LoadModel(Kratos.Parameters("""{
            "checkpoint_file" : "%s",
            "checkpoint_type" : "torchscript",
            "device"          : "cpu"
        }""" % self.checkpoint))
        second = training_utils.TrainModel(restored, self.dataset, Kratos.Parameters("""{
            "epochs" : 5, "device" : "cpu", "seed" : 0,
            "warm_restart" : { "shrink" : 0.8, "perturb" : 0.05 }
        }"""))
        self.assertEqual(len(second), 5)
        self.assertTrue(all(numpy.isfinite(value) for value in second))
        self.assertLess(min(second), first[0])


if __name__ == '__main__':
    KratosUnittest.main()
