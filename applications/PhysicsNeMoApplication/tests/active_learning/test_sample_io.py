from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.active_learning import sample_io
from KratosMultiphysics.PhysicsNeMoApplication.processes.export.dataset_export_process import DatasetExportProcess


class TestApplyParameterOverrides(KratosUnittest.TestCase):
    def setUp(self):
        self.parameters = Kratos.Parameters("""{
            "solver_settings": {
                "time_step": 0.1,
                "max_iterations": 10,
                "scheme": "backward_euler",
                "use_line_search": false
            },
            "processes": [
                { "Parameters": { "modulus": 1.0 } },
                { "Parameters": { "modulus": 2.0 } }
            ]
        }""")

    def test_ScalarOverrides(self):
        sample_io.ApplyParameterOverrides(self.parameters, {
            "solver_settings/time_step": 0.5,
            "solver_settings/max_iterations": 25,
            "solver_settings/scheme": "bdf2",
            "solver_settings/use_line_search": True,
        })
        self.assertAlmostEqual(self.parameters["solver_settings"]["time_step"].GetDouble(), 0.5)
        self.assertEqual(self.parameters["solver_settings"]["max_iterations"].GetInt(), 25)
        self.assertEqual(self.parameters["solver_settings"]["scheme"].GetString(), "bdf2")
        self.assertTrue(self.parameters["solver_settings"]["use_line_search"].GetBool())

    def test_ArrayIndexOverride(self):
        sample_io.ApplyParameterOverrides(self.parameters, {"processes/1/Parameters/modulus": 42.0})
        self.assertAlmostEqual(self.parameters["processes"][1]["Parameters"]["modulus"].GetDouble(), 42.0)
        self.assertAlmostEqual(self.parameters["processes"][0]["Parameters"]["modulus"].GetDouble(), 1.0)

    def test_IntIntoDoubleIsAccepted(self):
        sample_io.ApplyParameterOverrides(self.parameters, {"solver_settings/time_step": 1})
        self.assertAlmostEqual(self.parameters["solver_settings"]["time_step"].GetDouble(), 1.0)

    def test_MissingPathRaises(self):
        with self.assertRaisesRegex(ValueError, "does not exist"):
            sample_io.ApplyParameterOverrides(self.parameters, {"solver_settings/typo": 1.0})

    def test_TypeMismatchRaises(self):
        with self.assertRaisesRegex(ValueError, "Type mismatch"):
            sample_io.ApplyParameterOverrides(self.parameters, {"solver_settings/scheme": 3.0})

    def test_OutOfRangeIndexRaises(self):
        with self.assertRaisesRegex(ValueError, "out of range"):
            sample_io.ApplyParameterOverrides(self.parameters, {"processes/7/Parameters/modulus": 1.0})


class TestLoadFieldsFromNpzDirectory(KratosUnittest.TestCase):
    def setUp(self):
        self.output_path = Path("test_sample_io_dataset")
        self.model = Kratos.Model()
        model_part = self.model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        for i in range(4):
            node = model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.PRESSURE, 10.0 * (i + 1))

        settings = Kratos.Parameters("""{
            "model_part_name" : "Main",
            "list_of_fields"  : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
            "output_path"     : "test_sample_io_dataset",
            "output_interval" : 1
        }""")
        self.export = DatasetExportProcess(self.model, settings)
        self.export.ExecuteInitialize()
        for step in (1, 2):
            model_part.ProcessInfo[Kratos.STEP] = step
            model_part.ProcessInfo[Kratos.TIME] = 0.5 * step
            self.export.ExecuteFinalizeSolutionStep()

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))

    def test_LastStepOnly(self):
        fields, metadata = sample_io.LoadFieldsFromNpzDirectory(self.output_path)
        self.assertEqual(list(fields), ["PRESSURE__node_historical"])
        self.assertTrue(numpy.allclose(fields["PRESSURE__node_historical"], [10.0, 20.0, 30.0, 40.0]))
        self.assertEqual(int(metadata["STEP"][0]), 2)

    def test_AllSteps(self):
        fields, metadata = sample_io.LoadFieldsFromNpzDirectory(self.output_path, last_step_only=False)
        self.assertEqual(fields["PRESSURE__node_historical"].shape, (2, 4))
        self.assertEqual(metadata["STEP"].tolist(), [1, 2])

    def test_EmptyDirectoryRaises(self):
        empty = self.output_path / "empty"
        empty.mkdir()
        with self.assertRaises(FileNotFoundError):
            sample_io.LoadFieldsFromNpzDirectory(empty)


if __name__ == '__main__':
    KratosUnittest.main()
