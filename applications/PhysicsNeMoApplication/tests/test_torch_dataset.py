from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.training import torch_dataset
from KratosMultiphysics.PhysicsNeMoApplication.processes.export.dataset_export_process import DatasetExportProcess

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestTorchDataset(KratosUnittest.TestCase):
    def setUp(self):
        self.output_path = Path("test_torch_dataset")
        self.model = Kratos.Model()
        model_part = self.model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        for i in range(5):
            node = model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.PRESSURE, i + 1.0)
            node.SetSolutionStepValue(Kratos.VELOCITY, [i + 1.0, 0.0, 2.0])

        settings = Kratos.Parameters("""{
            "model_part_name" : "Main",
            "list_of_fields"  : [
                { "variable_name" : "VELOCITY", "data_location" : "node_historical" },
                { "variable_name" : "PRESSURE", "data_location" : "node_historical" }
            ],
            "output_path"     : "test_torch_dataset",
            "output_interval" : 1
        }""")
        export = DatasetExportProcess(self.model, settings)
        export.ExecuteInitialize()
        for step in (1, 2, 3):
            model_part.ProcessInfo[Kratos.STEP] = step
            model_part.ProcessInfo[Kratos.TIME] = float(step)
            export.ExecuteFinalizeSolutionStep()

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))

    def test_DatasetItems(self):
        dataset = torch_dataset.CreateNpzDataset(
            self.output_path,
            input_keys=["VELOCITY__node_historical"],
            output_keys=["PRESSURE__node_historical"])

        self.assertEqual(len(dataset), 3)
        inputs, outputs = dataset[0]
        self.assertEqual(list(inputs.shape), [5, 3])
        self.assertEqual(list(outputs.shape), [5, 1])
        self.assertEqual(inputs.dtype, torch.float32)
        self.assertTrue(numpy.allclose(outputs.numpy().ravel(), numpy.arange(1.0, 6.0)))

    def test_ConcatenatedInputKeys(self):
        dataset = torch_dataset.CreateNpzDataset(
            self.output_path,
            input_keys=["VELOCITY__node_historical", "PRESSURE__node_historical"],
            output_keys=["PRESSURE__node_historical"])
        inputs, _ = dataset[1]
        self.assertEqual(list(inputs.shape), [5, 4])

    def test_MissingKeyRaises(self):
        dataset = torch_dataset.CreateNpzDataset(
            self.output_path,
            input_keys=["TEMPERATURE__node_historical"],
            output_keys=["PRESSURE__node_historical"])
        with self.assertRaises(KeyError):
            dataset[0]

    def test_EmptyDirectoryRaises(self):
        empty = self.output_path / "empty"
        empty.mkdir()
        with self.assertRaises(FileNotFoundError):
            torch_dataset.CreateNpzDataset(empty, [], [])


if __name__ == '__main__':
    KratosUnittest.main()
