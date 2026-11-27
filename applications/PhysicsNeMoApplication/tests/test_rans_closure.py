"""A learned turbulence closure on a RANS channel flow.

RANSApplication is compiled in the reference build and nothing here had
driven it. The quantity a closure surrogate predicts is the modelled
TURBULENT_VISCOSITY: the two transport equations produce it, and a
surrogate learns it from the resolved flow instead.
"""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))

have_rans = kratos_utils.CheckIfApplicationsAvailable(
    "RANSApplication", "FluidDynamicsApplication", "LinearSolversApplication")

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False


@KratosUnittest.skipUnless(have_rans, "Missing RANS/FluidDynamics/LinearSolvers.")
class TestRansCase(KratosUnittest.TestCase):
    def setUp(self):
        import rans_case

        if not rans_case.IsAvailable():
            self.skipTest("RANSApplication channel-flow data is not present.")

    def test_TheChannelSolvesWithItsClosure(self):
        import rans_case

        features, viscosity, model_part = rans_case.SolveRans(Kratos.Model())
        self.assertEqual(features.shape, (model_part.NumberOfNodes(), 4))
        self.assertTrue(numpy.isfinite(features).all())
        self.assertTrue(numpy.isfinite(viscosity).all())
        self.assertGreater(numpy.abs(viscosity).max(), 0.0)
        # a closure that is constant across the channel is not a closure
        self.assertGreater(float(viscosity.std()), 0.0)

    def test_TheTurbulenceVariablesComeFromTheApplication(self):
        """Worth pinning: TURBULENT_KINETIC_ENERGY and its dissipation rate
        are registered by RANSApplication, not by the core, so reaching them
        through the Kratos namespace raises rather than returning a value."""
        import rans_case  # noqa: F401  (imports the application)

        with self.assertRaises(AttributeError):
            Kratos.TURBULENT_KINETIC_ENERGY
        self.assertIsNotNone(
            Kratos.KratosGlobals.GetVariable("TURBULENT_KINETIC_ENERGY"))

    def test_TheParametersAreATemplate(self):
        """The owning application ships placeholders its own driver fills
        in per test, so this case substitutes them rather than shipping a
        second copy of the file."""
        import rans_case

        text = rans_case._PARAMETERS_FILE.read_text()
        self.assertIn("<STABILIZATION_METHOD>", text)
        for placeholder in rans_case._SUBSTITUTIONS:
            self.assertIn(placeholder, text)


@KratosUnittest.skipUnless(have_rans and have_torch,
                           "Missing RANS applications or torch.")
class TestLearnedClosure(KratosUnittest.TestCase):
    def setUp(self):
        import rans_case

        if not rans_case.IsAvailable():
            self.skipTest("RANSApplication channel-flow data is not present.")
        self.features, self.viscosity, self.model_part = rans_case.SolveRans(Kratos.Model())

    def test_TheClosureIsLearnableFromTheResolvedFlow(self):
        """The modelled viscosity is a FUNCTION of the resolved quantities,
        which is what makes replacing the transport equations conceivable at
        all. Held out across nodes, the surrogate must beat predicting the
        channel's mean viscosity."""
        from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils

        features = torch.tensor(self.features, dtype=torch.float32)
        target = torch.tensor(self.viscosity, dtype=torch.float32).reshape(-1, 1)
        feature_scale = features.abs().amax(dim=0).clamp(min=1e-12)
        target_scale = target.abs().max().clamp(min=1e-12)
        features, target = features / feature_scale, target / target_scale

        generator = torch.Generator().manual_seed(0)
        order = torch.randperm(features.shape[0], generator=generator)
        split = int(0.7 * features.shape[0])
        train_rows, test_rows = order[:split], order[split:]

        torch.manual_seed(0)
        model = torch.nn.Sequential(
            torch.nn.Linear(4, 32), torch.nn.Tanh(), torch.nn.Linear(32, 1))
        history = training_utils.TrainModel(
            model,
            torch.utils.data.TensorDataset(features[train_rows], target[train_rows]),
            Kratos.Parameters("""{
                "epochs"        : 300,
                "batch_size"    : 16,
                "learning_rate" : 1e-2,
                "device"        : "cpu",
                "seed"          : 0
            }"""))
        self.assertLess(history[-1], history[0])

        with torch.no_grad():
            predicted = model(features[test_rows])
        surrogate_error = float((predicted - target[test_rows]).square().mean())
        mean_error = float((target[train_rows].mean() - target[test_rows]).square().mean())
        Kratos.Logger.PrintInfo(
            "TestLearnedClosure",
            f"held-out MSE {surrogate_error:.3e} against the mean predictor's {mean_error:.3e}")
        self.assertLess(surrogate_error, mean_error)

    def test_TheClosureFieldDeploysThroughTheOrdinaryProcess(self):
        """Nothing turbulence-specific is needed to write a closure back:
        TURBULENT_VISCOSITY is an ordinary nodal field to the bridge."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
        from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import (
            GetTensorAdaptor)

        adaptor = GetTensorAdaptor(
            self.model_part, "node_historical",
            Kratos.KratosGlobals.GetVariable("TURBULENT_VISCOSITY"))
        gathered = torch_bridge.KratosTensorToTorch(adaptor)
        numpy.testing.assert_allclose(
            gathered.numpy().reshape(-1), self.viscosity, atol=1e-12)


if __name__ == '__main__':
    KratosUnittest.main()
