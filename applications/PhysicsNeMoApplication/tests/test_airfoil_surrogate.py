"""An aerodynamic surrogate on compressible potential flow.

CompressiblePotentialFlowApplication is compiled and `adjoint_bridge`
already dispatches its response functions, but no case had ever been run
behind that dispatch. This sweeps the flight condition of the application's
small NACA 0012 case and learns the field it produces.
"""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))

have_aero = kratos_utils.CheckIfApplicationsAvailable(
    "CompressiblePotentialFlowApplication", "LinearSolversApplication")

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

_MACH_NUMBERS = (0.30, 0.40, 0.50, 0.60, 0.70)


@KratosUnittest.skipUnless(have_aero, "Missing CompressiblePotentialFlow/LinearSolvers.")
class TestAirfoilCase(KratosUnittest.TestCase):
    def setUp(self):
        import airfoil_case

        if not airfoil_case.IsAvailable():
            self.skipTest("CompressiblePotentialFlow test data is not present.")

    def test_TheFlightConditionDrivesTheSolution(self):
        import airfoil_case

        low, model_part = airfoil_case.SolveAirfoil(Kratos.Model(), mach_infinity=0.30)
        high, _ = airfoil_case.SolveAirfoil(Kratos.Model(), mach_infinity=0.70)
        self.assertTrue(numpy.isfinite(low).all())
        self.assertGreater(numpy.abs(high - low).max(), 1e-9)
        # compressibility stiffens the flow: a higher Mach number raises the
        # potential's magnitude on this case
        self.assertGreater(numpy.abs(high).max(), numpy.abs(low).max())

    def test_TheAirfoilSkinIsItsOwnModelPart(self):
        """Where a point-cloud surrogate would live."""
        import airfoil_case

        model = Kratos.Model()
        airfoil_case.SolveAirfoil(model, mach_infinity=0.5)
        body = airfoil_case.GetBodyModelPart(model)
        self.assertGreater(body.NumberOfNodes(), 0)
        self.assertLess(body.NumberOfNodes(),
                        airfoil_case.GetModelPart(model).NumberOfNodes())

    def test_ThisCaseHasNoAngleOfAttackToSet(self):
        """Worth pinning because it is a silent trap: the compressible NACA
        case is parameterized by Mach number alone, so a caller sweeping an
        angle of attack would get identical solves and conclude the
        surrogate had learned something."""
        import airfoil_case

        with self.assertRaisesRegex(ValueError, "angle_of_attack"):
            airfoil_case.CreateAirfoilAnalysis(
                Kratos.Model(), mach_infinity=0.5, angle_of_attack=5.0)


@KratosUnittest.skipUnless(have_aero and have_torch,
                           "Missing CompressiblePotentialFlow/LinearSolvers or torch.")
class TestAirfoilSurrogate(KratosUnittest.TestCase):
    """Mach number in, the whole potential field out - the operator shape
    `RomSurrogateProcess` covers through a basis, learned directly."""

    def setUp(self):
        import airfoil_case

        if not airfoil_case.IsAvailable():
            self.skipTest("CompressiblePotentialFlow test data is not present.")
        self.fields = []
        for mach in _MACH_NUMBERS:
            potentials, _ = airfoil_case.SolveAirfoil(Kratos.Model(), mach_infinity=mach)
            self.fields.append(potentials)
        self.fields = numpy.stack(self.fields)

    def test_TheSweepIsLearnable(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils

        machs = torch.tensor(_MACH_NUMBERS, dtype=torch.float32).reshape(-1, 1)
        fields = torch.tensor(self.fields, dtype=torch.float32)
        scale = fields.abs().max()
        targets = fields / scale

        train = torch.utils.data.TensorDataset(machs[:-1], targets[:-1])
        torch.manual_seed(0)
        model = torch.nn.Sequential(
            torch.nn.Linear(1, 32), torch.nn.Tanh(),
            torch.nn.Linear(32, targets.shape[1]))
        history = training_utils.TrainModel(model, train, Kratos.Parameters("""{
            "epochs"        : 400,
            "batch_size"    : 4,
            "learning_rate" : 1e-2,
            "device"        : "cpu",
            "shuffle"       : false,
            "seed"          : 0
        }"""))
        self.assertLess(history[-1], history[0])

        # the held-out flight condition, against predicting the sweep's mean
        with torch.no_grad():
            predicted = model(machs[-1:])
        surrogate_error = float((predicted - targets[-1:]).square().mean())
        mean_error = float((targets[:-1].mean(dim=0, keepdim=True)
                            - targets[-1:]).square().mean())
        Kratos.Logger.PrintInfo(
            "TestAirfoilSurrogate",
            f"held-out MSE {surrogate_error:.3e} against the mean predictor's {mean_error:.3e}")
        self.assertLess(surrogate_error, mean_error)


@KratosUnittest.skipUnless(have_aero, "Missing CompressiblePotentialFlow/LinearSolvers.")
class TestAirfoilAdjointDispatch(KratosUnittest.TestCase):
    """The adjoint bridge's aero entry, on a compiled application.

    The dispatch table has always named this application; what was missing
    was any evidence it resolves. A full adjoint cross-validation needs the
    separate adjoint solver and primal-restart files the application ships
    for its own regression, which this case does not set up - see the
    roadmap."""

    def test_TheResponseFactoryResolves(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import adjoint_bridge

        application, module_path = adjoint_bridge._RESPONSE_FACTORIES[
            "compressible_potential_flow"]
        self.assertIn("CompressiblePotentialFlowApplication", application)
        import importlib

        factory = importlib.import_module(module_path)
        self.assertTrue(hasattr(factory, "CreateResponseFunction"))

    def test_AnUnknownResponseTypeIsReported(self):
        import importlib

        from KratosMultiphysics.PhysicsNeMoApplication.bridges import adjoint_bridge

        _, module_path = adjoint_bridge._RESPONSE_FACTORIES["compressible_potential_flow"]
        factory = importlib.import_module(module_path)
        settings = Kratos.Parameters('{ "response_type" : "drag" }')
        with self.assertRaises(Exception):
            factory.CreateResponseFunction("r", settings, Kratos.Model())


if __name__ == '__main__':
    KratosUnittest.main()
