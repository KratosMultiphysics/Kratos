from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
from KratosMultiphysics.PhysicsNeMoApplication.deployment import uncertainty_utils
from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.nn import ConcreteDropout  # noqa: F401
    have_concrete_dropout = True
except ImportError:
    have_concrete_dropout = False


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestMonteCarloPredict(KratosUnittest.TestCase):
    def test_MeanAndStdOverDropoutSamples(self):
        torch.manual_seed(0)
        model = torch.nn.Sequential(
            torch.nn.Linear(2, 16), torch.nn.Dropout(p=0.5), torch.nn.Linear(16, 1)).double()
        model.eval()
        x = torch.rand(10, 2, dtype=torch.float64)

        mean, std = uncertainty_utils.MonteCarloPredict(
            model, lambda m: m(x), num_samples=16, seed=0)
        self.assertEqual(list(mean.shape), [10, 1])
        self.assertEqual(list(std.shape), [10, 1])
        self.assertGreater(float(std.abs().max()), 0.0)
        # the sampling must not leave the model in train mode
        self.assertFalse(model.training)
        for module in model.modules():
            self.assertFalse(module.training)

    def test_SeededSamplingIsReproducible(self):
        torch.manual_seed(0)
        model = torch.nn.Sequential(
            torch.nn.Linear(2, 8), torch.nn.Dropout(p=0.3), torch.nn.Linear(8, 1)).double()
        x = torch.rand(5, 2, dtype=torch.float64)
        mean_a, std_a = uncertainty_utils.MonteCarloPredict(model, lambda m: m(x), 8, seed=42)
        mean_b, std_b = uncertainty_utils.MonteCarloPredict(model, lambda m: m(x), 8, seed=42)
        self.assertTrue(torch.equal(mean_a, mean_b))
        self.assertTrue(torch.equal(std_a, std_b))

    def test_ModelWithoutDropoutRaises(self):
        model = torch.nn.Linear(2, 1).double()
        with self.assertRaisesRegex(ValueError, "dropout-like"):
            uncertainty_utils.MonteCarloPredict(model, lambda m: m(torch.zeros(3, 2)), 4)

    def test_TooFewSamplesRaises(self):
        model = torch.nn.Sequential(torch.nn.Dropout(0.5))
        with self.assertRaisesRegex(ValueError, "num_samples"):
            uncertainty_utils.MonteCarloPredict(model, lambda m: m(torch.zeros(3, 2)), 1)


@KratosUnittest.skipUnless(have_torch and have_concrete_dropout,
                           "Missing required python modules: torch, physicsnemo.")
class TestConcreteDropoutTraining(KratosUnittest.TestCase):
    def test_TrainModelWithRegularization(self):
        from physicsnemo.nn import ConcreteDropout

        torch.manual_seed(0)
        model = torch.nn.Sequential(
            torch.nn.Linear(2, 16), ConcreteDropout(in_features=16), torch.nn.Linear(16, 1))
        inputs = torch.rand(64, 2)
        dataset = torch.utils.data.TensorDataset(inputs, inputs.sum(dim=1, keepdim=True))

        history = training_utils.TrainModel(model, dataset, Kratos.Parameters("""{
            "epochs"                      : 3,
            "seed"                        : 0,
            "device"                      : "cpu",
            "concrete_dropout_reg_weight" : 1e-4
        }"""))
        self.assertEqual(len(history), 3)

        rates = uncertainty_utils.GetConcreteDropoutRates(model)
        self.assertEqual(len(rates), 1)
        for rate in rates.values():
            self.assertGreater(rate, 0.0)
            self.assertLess(rate, 1.0)
        # MC inference works through the ConcreteDropout layer
        mean, std = uncertainty_utils.MonteCarloPredict(
            model, lambda m: m(inputs[:5]), 8, seed=0)
        self.assertGreater(float(std.abs().max()), 0.0)


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestProcessUncertainty(KratosUnittest.TestCase):
    """The "uncertainty" block of InferenceProcess: mc_dropout and ensemble."""

    def setUp(self):
        self.checkpoints = []
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(
            self.model, "Main", divisions=2,
            historical_variables=(Kratos.PRESSURE, Kratos.TEMPERATURE, Kratos.NODAL_PAUX))
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0 + node.X)

    def tearDown(self):
        for path in self.checkpoints:
            KratosUtilities.DeleteFileIfExisting(str(path))

    def _SaveScale(self, path, scale):
        class Scale(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.scale = scale

            def forward(self, x):
                return self.scale * x

        torch.jit.script(Scale()).save(str(path))
        self.checkpoints.append(path)

    def test_EnsembleMeanAndStd(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
        self._SaveScale(Path("test_unc_ensemble_a.pt"), 2.0)
        self._SaveScale(Path("test_unc_ensemble_b.pt"), 4.0)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {
                    "checkpoint_files" : ["test_unc_ensemble_a.pt", "test_unc_ensemble_b.pt"],
                    "device"           : "cpu"
                },
                "uncertainty"     : {
                    "method"             : "ensemble",
                    "uncertainty_fields" : [ { "variable_name" : "NODAL_PAUX", "data_location" : "node_historical" } ]
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""")
        process = inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        for node in self.model_part.Nodes:
            pressure = 1.0 + node.X
            # members predict 2p and 4p: mean 3p, std |p| * sqrt(2)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.TEMPERATURE), 3.0 * pressure, places=10)
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.NODAL_PAUX),
                pressure * (2.0 ** 0.5), places=10)

    def test_EnsembleWithoutFilesRaises(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
        self._SaveScale(Path("test_unc_single.pt"), 2.0)
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : { "checkpoint_file" : "test_unc_single.pt", "device" : "cpu" },
                "uncertainty"     : { "method" : "ensemble" },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""")
        process = inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        with self.assertRaisesRegex(ValueError, "checkpoint_files"):
            process.ExecuteFinalizeSolutionStep()

    def test_McDropoutThroughProcess(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
        torch.manual_seed(0)
        dropout_model = torch.nn.Sequential(
            torch.nn.Linear(1, 16), torch.nn.Dropout(p=0.5), torch.nn.Linear(16, 1)).double()
        checkpoint = Path("test_unc_mc_model.pt")
        torch.jit.script(dropout_model).save(str(checkpoint))
        self.checkpoints.append(checkpoint)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : { "checkpoint_file" : "test_unc_mc_model.pt", "device" : "cpu" },
                "uncertainty"     : {
                    "method"             : "mc_dropout",
                    "num_samples"        : 8,
                    "seed"               : 0,
                    "uncertainty_fields" : [ { "variable_name" : "NODAL_PAUX", "data_location" : "node_historical" } ]
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""")
        process = inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        stds = numpy.array([
            node.GetSolutionStepValue(Kratos.NODAL_PAUX) for node in self.model_part.Nodes])
        self.assertTrue(numpy.isfinite(stds).all())
        self.assertGreater(numpy.abs(stds).max(), 0.0)

    def test_UnknownMethodRaises(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {},
                "uncertainty"     : { "method" : "bootstrap" },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""")
        with self.assertRaisesRegex(ValueError, "uncertainty method"):
            inference_process.Factory(settings, self.model)


if __name__ == '__main__':
    KratosUnittest.main()
