"""Integration tests driven by REAL Kratos solves.

Conditioned on the compiled applications following the
CoSimulationApplication pattern: module-level availability flags from
kratos_utilities.CheckIfApplicationsAvailable + skipUnless decorators. Every
test here runs an actual ConvectionDiffusionApplication stationary solve.
"""

import os
import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

have_convection_diffusion = kratos_utils.CheckIfApplicationsAvailable(
    "ConvectionDiffusionApplication", "LinearSolversApplication")

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import physicsnemo.metrics.general.mse
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))


def _RunRealSolve(model, conductivity, divisions=8):
    import thermal_case
    analysis = thermal_case.CreateThermalAnalysis(
        model, conductivity=conductivity, heat_flux=1.0, divisions=divisions)
    analysis.Run()
    return model["ThermalModelPart"]


@KratosUnittest.skipUnless(have_convection_diffusion,
                           "Missing required applications: ConvectionDiffusionApplication, LinearSolversApplication.")
class TestRealSolverDatasetExport(KratosUnittest.TestCase):
    def setUp(self):
        self.output_path = Path("test_real_solver_dataset")

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting(str(self.output_path))

    def test_ExportedFieldIsTheSolversSolution(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.export.dataset_export_process import DatasetExportProcess

        model = Kratos.Model()
        model_part = _RunRealSolve(model, conductivity=2.0)
        model_part.ProcessInfo[Kratos.STEP] = 1

        settings = Kratos.Parameters("""{
            "model_part_name" : "ThermalModelPart",
            "list_of_fields"  : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
            "output_path"     : "test_real_solver_dataset",
            "output_interval" : 1
        }""")
        export = DatasetExportProcess(model, settings)
        export.ExecuteInitialize()
        export.ExecuteFinalizeSolutionStep()

        with numpy.load(self.output_path / "sample_1.npz") as data:
            temperature = data["TEMPERATURE__node_historical"]
        # boundary held at zero, interior heated: physical, non-trivial field
        self.assertEqual(temperature.shape[0], model_part.NumberOfNodes())
        self.assertTrue(numpy.isfinite(temperature).all())
        self.assertGreater(temperature.max(), 0.01)
        self.assertAlmostEqual(temperature.min(), 0.0, places=10)
        # Poisson max on the unit square: ~0.0737 * f / k
        self.assertAlmostEqual(temperature.max(), 0.0737 / 2.0, delta=0.005)


@KratosUnittest.skipUnless(have_convection_diffusion,
                           "Missing required applications: ConvectionDiffusionApplication, LinearSolversApplication.")
class TestRealSolverActiveLearning(KratosUnittest.TestCase):
    """SubprocessBackend labeling with an actual external Kratos solve."""

    def setUp(self):
        self.working_directory = Path("test_real_solver_al_cases")
        self._previous_pythonpath = os.environ.get("PYTHONPATH", "")
        os.environ["PYTHONPATH"] = str(_CASES_DIR) + os.pathsep + self._previous_pythonpath

    def tearDown(self):
        os.environ["PYTHONPATH"] = self._previous_pythonpath
        kratos_utils.DeleteDirectoryIfExisting(str(self.working_directory))

    def test_LabelingWithRealSolves(self):
        from KratosMultiphysics.PhysicsNeMoApplication.active_learning.sample_io import KratosALSample
        from KratosMultiphysics.PhysicsNeMoApplication.active_learning.execution_backends.subprocess_backend import SubprocessBackend

        settings = Kratos.Parameters("""{
            "template_directory" : "",
            "working_directory"  : "test_real_solver_al_cases",
            "timeout_seconds"    : 300,
            "max_retries"        : 0
        }""")
        settings["template_directory"].SetString(str(_CASES_DIR / "thermal_template"))
        settings.AddEmptyArray("run_command")
        settings["run_command"].Append(sys.executable)
        settings["run_command"].Append("MainKratos.py")
        backend = SubprocessBackend(settings)

        soft = backend.RunCase(KratosALSample("soft", parameters={"case_settings/conductivity": 1.0}))
        stiff = backend.RunCase(KratosALSample("stiff", parameters={"case_settings/conductivity": 4.0}))

        soft_temperature = soft.fields["TEMPERATURE__node_historical"]
        stiff_temperature = stiff.fields["TEMPERATURE__node_historical"]
        # Quadrupled conductivity -> roughly quartered temperature (exactly 4
        # for the pure Poisson problem; the stabilized element makes the
        # discrete scaling only approximately linear in k).
        self.assertAlmostEqual(soft_temperature.max() / stiff_temperature.max(), 4.0, delta=0.5)


@KratosUnittest.skipUnless(have_convection_diffusion and have_torch and have_physicsnemo,
                           "Missing ConvectionDiffusionApplication/LinearSolversApplication, torch or physicsnemo.")
class TestRealSolverSurrogateValidation(KratosUnittest.TestCase):
    """Train on real solves, deploy, and validate with ValidationMetricsProcess."""

    def setUp(self):
        self.checkpoint = Path("test_real_solver_surrogate.pt")
        self.report = Path("test_real_solver_metrics.json")

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting(str(self.checkpoint))
        kratos_utils.DeleteFileIfExisting(str(self.report))

    def test_SurrogateOfRealSolves(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
        from KratosMultiphysics.PhysicsNeMoApplication.processes import validation_metrics_process
        # Training sweep: real solves at several conductivities. Node-local
        # inputs (x, y, k) -> solved TEMPERATURE.
        train_inputs, train_targets = [], []
        for conductivity in (0.5, 1.0, 2.0, 4.0):
            model = Kratos.Model()
            model_part = _RunRealSolve(model, conductivity)
            train_inputs.append(numpy.array(
                [[node.X, node.Y, conductivity] for node in model_part.Nodes]))
            train_targets.append(numpy.array(
                [[node.GetSolutionStepValue(Kratos.TEMPERATURE)] for node in model_part.Nodes]))
        inputs = torch.from_numpy(numpy.concatenate(train_inputs)).double()
        targets = torch.from_numpy(numpy.concatenate(train_targets)).double()

        surrogate = torch.nn.Sequential(
            torch.nn.Linear(3, 32), torch.nn.Tanh(),
            torch.nn.Linear(32, 32), torch.nn.Tanh(),
            torch.nn.Linear(32, 1)).double()
        optimizer = torch.optim.Adam(surrogate.parameters(), lr=5e-3)
        for epoch in range(500):
            optimizer.zero_grad()
            loss = torch.nn.functional.mse_loss(surrogate(inputs), targets)
            loss.backward()
            optimizer.step()
        torch.jit.script(surrogate).save(str(self.checkpoint))

        # Held-out real solve at an unseen conductivity. Surrogate inputs are
        # packed into non-historical NODAL_VAUX (x, y, k); the prediction
        # lands in non-historical NODAL_PAUX.
        model = Kratos.Model()
        model_part = _RunRealSolve(model, conductivity=1.5)
        for node in model_part.Nodes:
            node.SetValue(Kratos.NODAL_VAUX, [node.X, node.Y, 1.5])
            node.SetValue(Kratos.NODAL_PAUX, 0.0)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "ThermalModelPart",
                "model_settings"  : { "checkpoint_file" : "test_real_solver_surrogate.pt", "device" : "cpu" },
                "input_fields"    : [ { "variable_name" : "NODAL_VAUX", "data_location" : "node_non_historical" } ],
                "output_fields"   : [ { "variable_name" : "NODAL_PAUX", "data_location" : "node_non_historical" } ]
            }
        }""")
        process = inference_process.Factory(settings, model)
        model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        metrics_settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"     : "ThermalModelPart",
                "list_of_comparisons" : [{
                    "predicted_variable" : "NODAL_PAUX",
                    "predicted_location" : "node_non_historical",
                    "reference_variable" : "TEMPERATURE",
                    "reference_location" : "node_historical",
                    "metrics"            : ["rmse", "max_abs_error"]
                }],
                "output_file"         : "test_real_solver_metrics.json"
            }
        }""")
        metrics = validation_metrics_process.Factory(metrics_settings, model)
        metrics.ExecuteFinalizeSolutionStep()
        metrics.ExecuteFinalize()

        values = metrics.history[0]["NODAL_PAUX_vs_TEMPERATURE"]
        # Loose threshold: the surrogate must beat the trivial zero predictor
        # (field max ~0.05) by a clear margin on the held-out real solve.
        self.assertLess(values["rmse"], 0.015)


if __name__ == '__main__':
    KratosUnittest.main()
