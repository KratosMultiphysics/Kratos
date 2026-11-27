"""Integration tests driven by REAL StructuralMechanics solves.

The structural counterpart of test_real_solver_integration.py: a
crash/deformation-surrogate pattern - load parameter in, displacement field
out - trained and deployed on actual linear-static cantilever solves.
Availability-gated on the compiled applications.
"""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

have_structural = kratos_utils.CheckIfApplicationsAvailable(
    "StructuralMechanicsApplication", "LinearSolversApplication")

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))


def _RunRealSolve(model, tip_load, divisions=6):
    import structural_case
    analysis = structural_case.CreateStructuralAnalysis(
        model, tip_load=tip_load, divisions=divisions)
    analysis.Run()
    return model["StructuralModelPart"]


@KratosUnittest.skipUnless(have_structural,
                           "Missing required applications: StructuralMechanicsApplication, LinearSolversApplication.")
class TestStructuralSolverCase(KratosUnittest.TestCase):
    def test_CantileverDeflectsDownwardAndLinearly(self):
        import structural_case

        model_part_1 = _RunRealSolve(Kratos.Model(), tip_load=1.0e6)
        deflection_1 = structural_case.GetTipDeflection(model_part_1)
        self.assertLess(deflection_1, 0.0)  # downward load, downward deflection

        model_part_2 = _RunRealSolve(Kratos.Model(), tip_load=2.0e6)
        deflection_2 = structural_case.GetTipDeflection(model_part_2)
        # linear problem: the response is exactly proportional to the load
        self.assertAlmostEqual(deflection_2 / deflection_1, 2.0, places=8)

    def test_ClampedEdgeStaysFixed(self):
        model_part = _RunRealSolve(Kratos.Model(), tip_load=1.0e6)
        for node in model_part.Nodes:
            if abs(node.X0) < 1e-8:
                self.assertEqual(node.GetSolutionStepValue(Kratos.DISPLACEMENT_X), 0.0)
                self.assertEqual(node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y), 0.0)


@KratosUnittest.skipUnless(have_structural and have_torch,
                           "Missing required applications/modules: StructuralMechanicsApplication, "
                           "LinearSolversApplication, torch.")
class TestStructuralSolverSurrogate(KratosUnittest.TestCase):
    """Deformation surrogate: load parameter -> DISPLACEMENT field."""

    def setUp(self):
        self.checkpoint = Path("test_structural_surrogate.pt")

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting(str(self.checkpoint))

    @staticmethod
    def _GatherDisplacements(model_part):
        return numpy.array([
            list(node.GetSolutionStepValue(Kratos.DISPLACEMENT)) for node in model_part.Nodes])

    def test_SurrogateReproducesUnseenLoad(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
        # --- solve the training loads with the real solver ------------------
        train_loads = (0.5e6, 1.0e6, 2.0e6)
        samples = []
        for load in train_loads:
            model_part = _RunRealSolve(Kratos.Model(), tip_load=load)
            samples.append((load, self._GatherDisplacements(model_part)))
        n_nodes = samples[0][1].shape[0]

        # --- train the load -> displacement-field map -----------------------
        # The linear response is u(load) = load * u_unit; the surrogate learns
        # the per-node unit field u_unit from the solves via TrainModel.
        load_scale = 1.0e6

        class DeformationSurrogate(torch.nn.Module):
            def __init__(self, n_nodes):
                super().__init__()
                self.unit_displacements = torch.nn.Parameter(
                    torch.zeros((n_nodes, 3), dtype=torch.float64))

            def forward(self, load_fraction):  # (N, 1) -> (N, 3)
                return load_fraction * self.unit_displacements

        inputs = torch.stack([
            torch.full((n_nodes, 1), load / load_scale, dtype=torch.float64)
            for load, _ in samples])
        targets = torch.stack(
            [torch.tensor(displacements, dtype=torch.float64) for _, displacements in samples])
        dataset = torch.utils.data.TensorDataset(inputs, targets)

        surrogate = DeformationSurrogate(n_nodes)
        history = training_utils.TrainModel(surrogate, dataset, Kratos.Parameters("""{
            "epochs"        : 400,
            "batch_size"    : 3,
            "learning_rate" : 1e-3,
            "device"        : "cpu",
            "seed"          : 0
        }"""))
        self.assertLess(history[-1], history[0])
        torch.jit.script(surrogate).save(str(self.checkpoint))

        # --- deploy on a REAL solve of an unseen load -----------------------
        test_load = 1.5e6
        model = Kratos.Model()
        model_part = _RunRealSolve(model, tip_load=test_load)
        reference = self._GatherDisplacements(model_part)

        # carry the load parameter in a pre-added scalar nodal variable
        for node in model_part.Nodes:
            node.SetSolutionStepValue(Kratos.POSITIVE_FACE_PRESSURE, test_load / load_scale)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "StructuralModelPart",
                "model_settings"  : { "checkpoint_file" : "%s", "device" : "cpu" },
                "input_fields"    : [ { "variable_name" : "POSITIVE_FACE_PRESSURE", "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "MESH_DISPLACEMENT",      "data_location" : "node_non_historical" } ]
            }
        }""" % self.checkpoint)
        process = inference_process.Factory(settings, model)
        model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        predicted = numpy.array([
            list(node.GetValue(Kratos.MESH_DISPLACEMENT)) for node in model_part.Nodes])
        scale = numpy.abs(reference).max()
        self.assertGreater(scale, 0.0)
        # linear problem -> the fitted unit field makes the surrogate accurate
        self.assertLess(numpy.abs(predicted - reference).max() / scale, 1e-2)


if __name__ == '__main__':
    KratosUnittest.main()
