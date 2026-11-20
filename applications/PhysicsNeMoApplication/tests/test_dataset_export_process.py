from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.processes.export import dataset_export_process
class TestDatasetExportProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test")
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        for i in range(6):
            node = self.model_part.CreateNewNode(i + 1, i, 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.PRESSURE, i + 1.0)
            node.SetSolutionStepValue(Kratos.VELOCITY, [i + 1.0, 0.0, -1.0])
        self.output_path = Path("test_physics_nemo_dataset")

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))

    def test_DatasetExport(self):
        settings = Kratos.Parameters("""{
            "Parameters" : {
                "model_part_name" : "test",
                "list_of_fields"  : [
                    { "variable_name" : "PRESSURE", "data_location" : "node_historical" },
                    { "variable_name" : "VELOCITY", "data_location" : "node_historical" }
                ],
                "output_path"     : "test_physics_nemo_dataset",
                "output_interval" : 2
            }
        }""")
        process = dataset_export_process.Factory(settings, self.model)
        process.ExecuteInitialize()

        for step in range(1, 5):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            self.model_part.ProcessInfo[Kratos.TIME] = 0.1 * step
            process.ExecuteFinalizeSolutionStep()

        # interval = 2 over steps 1..4 -> samples at steps 2 and 4 only
        written = sorted(p.name for p in self.output_path.glob("*.npz"))
        self.assertEqual(written, ["sample_2.npz", "sample_4.npz"])

        with numpy.load(self.output_path / "sample_2.npz") as data:
            self.assertEqual(int(data["STEP"]), 2)
            self.assertAlmostEqual(float(data["TIME"]), 0.2)
            pressure = data["PRESSURE__node_historical"]
            velocity = data["VELOCITY__node_historical"]
            self.assertEqual(pressure.shape, (6,))
            self.assertEqual(velocity.shape, (6, 3))
            self.assertTrue(numpy.allclose(pressure, numpy.arange(1.0, 7.0)))
            self.assertTrue(numpy.allclose(velocity[:, 0], numpy.arange(1.0, 7.0)))


if __name__ == '__main__':
    KratosUnittest.main()
