from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import physicsnemo
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestTrainModel(KratosUnittest.TestCase):
    def _CreateDataset(self):
        # y = 3x - 1 with a fixed generator: exactly learnable by Linear(1, 1).
        generator = torch.Generator().manual_seed(0)
        x = torch.rand(256, 1, generator=generator, dtype=torch.float64)
        return torch.utils.data.TensorDataset(x, 3.0 * x - 1.0)

    def test_ConvergesOnLinearProblem(self):
        model = torch.nn.Linear(1, 1).double()
        history = training_utils.TrainModel(model, self._CreateDataset(), Kratos.Parameters("""{
            "epochs"        : 200,
            "batch_size"    : 64,
            "learning_rate" : 0.05,
            "device"        : "cpu",
            "seed"          : 0
        }"""))
        self.assertEqual(len(history), 200)
        self.assertLess(history[-1], 1e-3)
        self.assertLess(history[-1], history[0])
        self.assertAlmostEqual(float(model.weight), 3.0, places=2)
        self.assertAlmostEqual(float(model.bias), -1.0, places=2)

    def test_SgdAndL1Paths(self):
        model = torch.nn.Linear(1, 1).double()
        history = training_utils.TrainModel(model, self._CreateDataset(), Kratos.Parameters("""{
            "epochs"        : 5,
            "optimizer"     : "sgd",
            "loss"          : "l1",
            "device"        : "cpu",
            "seed"          : 0
        }"""))
        self.assertEqual(len(history), 5)

    def test_BadOptimizerRaises(self):
        with self.assertRaisesRegex(ValueError, "Unsupported optimizer"):
            training_utils.TrainModel(
                torch.nn.Linear(1, 1), self._CreateDataset(),
                Kratos.Parameters('{"optimizer": "lbfgs"}'))

    def test_BadLossRaises(self):
        with self.assertRaisesRegex(ValueError, "Unsupported loss"):
            training_utils.TrainModel(
                torch.nn.Linear(1, 1), self._CreateDataset(),
                Kratos.Parameters('{"loss": "huber"}'))


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestTrainModelCallbacks(KratosUnittest.TestCase):
    def test_CallbacksRunAfterEveryEpoch(self):
        generator = torch.Generator().manual_seed(0)
        x = torch.rand(32, 1, generator=generator, dtype=torch.float64)
        dataset = torch.utils.data.TensorDataset(x, 2.0 * x)
        model = torch.nn.Linear(1, 1).double()

        calls = []

        def monitor(epoch, callback_model, history):
            # invoked with the model in eval mode, inside no_grad: the
            # physics-informed residual-monitor contract
            self.assertFalse(callback_model.training)
            self.assertFalse(torch.is_grad_enabled())
            self.assertEqual(len(history), epoch + 1)
            calls.append(epoch)

        history = training_utils.TrainModel(model, dataset, Kratos.Parameters("""{
            "epochs"     : 3,
            "batch_size" : 16,
            "device"     : "cpu",
            "seed"       : 0
        }"""), epoch_callbacks=[monitor])
        self.assertEqual(calls, [0, 1, 2])
        self.assertEqual(len(history), 3)
        self.assertFalse(model.training)


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestSaveTrainedModel(KratosUnittest.TestCase):
    def setUp(self):
        self.checkpoint = Path("test_training_utils_model.pt")

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint) + ".card.json")

    def test_TorchScriptRoundTrip(self):
        model = torch.nn.Linear(2, 1).double()
        checkpoint_type = training_utils.SaveTrainedModel(
            model, self.checkpoint,
            card={"input_fields": [{"variable_name": "VELOCITY", "data_location": "node_historical"}]})
        self.assertEqual(checkpoint_type, "torchscript")

        settings = Kratos.Parameters("""{
            "checkpoint_file" : "test_training_utils_model.pt",
            "checkpoint_type" : "torchscript",
            "device"          : "cpu"
        }""")
        loaded, _ = model_registry.LoadModel(settings)
        self.assertEqual(list(loaded(torch.zeros(3, 2, dtype=torch.float64)).shape), [3, 1])

        card = model_registry.LoadModelCard(self.checkpoint)
        self.assertEqual(card["input_fields"][0]["variable_name"], "VELOCITY")

    @KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
    def test_PhysicsNemoModuleRoundTrip(self):
        from physicsnemo.models import FullyConnected
        model = FullyConnected(in_features=2, out_features=1, num_layers=1, layer_size=4)

        with self.assertRaisesRegex(ValueError, ".mdlus"):
            training_utils.SaveTrainedModel(model, "test_training_utils_model.pt")

        try:
            checkpoint_type = training_utils.SaveTrainedModel(model, "test_training_utils_model.mdlus")
            self.assertEqual(checkpoint_type, "physicsnemo")

            settings = Kratos.Parameters("""{
                "checkpoint_file" : "test_training_utils_model.mdlus",
                "checkpoint_type" : "physicsnemo",
                "device"          : "cpu"
            }""")
            loaded, _ = model_registry.LoadModel(settings)
            self.assertEqual(list(loaded(torch.zeros(3, 2)).shape), [3, 1])
        finally:
            KratosUtilities.DeleteFileIfExisting("test_training_utils_model.mdlus")


try:
    from torch.distributed.fsdp import fully_shard  # noqa: F401
    from torch.distributed.tensor import DTensor  # noqa: F401
    have_fsdp = have_torch and have_physicsnemo
except ImportError:
    have_fsdp = False


@KratosUnittest.skipUnless(have_fsdp,
                           "Missing torch FSDP2 / physicsnemo.")
class TestSaveShardedModel(KratosUnittest.TestCase):
    """SaveTrainedModel must not write DTensors.

    FSDP2 replaces parameters with DTensors; serializing those produces a
    checkpoint that reports success and then cannot be loaded. This
    reproduces at world_size=1, so no MPI is needed.
    """

    _PORT = "29671"

    def setUp(self):
        import os
        import torch.distributed as distributed
        self.previous = {key: os.environ.get(key)
                         for key in ("MASTER_ADDR", "MASTER_PORT", "RANK", "WORLD_SIZE")}
        os.environ.update(MASTER_ADDR="127.0.0.1", MASTER_PORT=self._PORT,
                          RANK="0", WORLD_SIZE="1")
        self.started_group = not distributed.is_initialized()
        if self.started_group:
            distributed.init_process_group("gloo", rank=0, world_size=1)
        self.checkpoint = Path("test_sharded_save.mdlus")

    def tearDown(self):
        import os
        import torch.distributed as distributed
        if self.started_group and distributed.is_initialized():
            distributed.destroy_process_group()
        for key, value in self.previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    @staticmethod
    def _Model():
        from physicsnemo.models.mlp.fully_connected import FullyConnected
        torch.manual_seed(0)
        return FullyConnected(in_features=4, out_features=2,
                              layer_size=16, num_layers=2)

    def test_ShardedCheckpointContainsNoDTensorsAndReloads(self):
        import io
        import zipfile
        from torch.distributed.fsdp import fully_shard
        from torch.distributed.tensor import DTensor
        import physicsnemo

        model = self._Model()
        reference = {name: value.detach().clone()
                     for name, value in model.state_dict().items()}
        fully_shard(model)
        self.assertTrue(any(isinstance(p, DTensor) for p in model.parameters()))

        checkpoint_type = training_utils.SaveTrainedModel(model, str(self.checkpoint))
        self.assertEqual(checkpoint_type, "physicsnemo")

        # what actually landed on disk: plain tensors, not DTensors
        with zipfile.ZipFile(str(self.checkpoint)) as archive:
            stored = torch.load(io.BytesIO(archive.read("model.pt")), weights_only=False)
        self.assertEqual({type(v).__name__ for v in stored.values()}, {"Tensor"})

        reloaded = physicsnemo.Module.from_checkpoint(str(self.checkpoint))
        for name, value in reloaded.state_dict().items():
            self.assertTrue(torch.allclose(reference[name], value),
                            msg=f"{name} differs from the pre-shard weights")

    def test_UnshardedModelsAreUnaffected(self):
        model = self._Model()
        reference = {name: value.detach().clone()
                     for name, value in model.state_dict().items()}
        self.assertEqual(
            training_utils.SaveTrainedModel(model, str(self.checkpoint)), "physicsnemo")
        import physicsnemo
        reloaded = physicsnemo.Module.from_checkpoint(str(self.checkpoint))
        for name, value in reloaded.state_dict().items():
            self.assertTrue(torch.allclose(reference[name], value))

    def test_GatherReshardsFirst(self):
        # a forward whose backward never runs leaves parameters unsharded as
        # plain Parameters; without a reshard the gather would see no DTensor
        # and silently fall through to the plain save path
        from torch.distributed.fsdp import fully_shard
        from torch.distributed.tensor import DTensor

        model = self._Model()
        fully_shard(model)
        model(torch.zeros(2, 4))                       # no backward
        self.assertFalse(any(isinstance(p, DTensor) for p in model.parameters()))

        gathered = training_utils._GatherShardedStateDict(model)
        self.assertIsNotNone(gathered)                 # reshard restored them
        self.assertEqual({type(v).__name__ for v in gathered.values()}, {"Tensor"})


if __name__ == '__main__':
    KratosUnittest.main()
