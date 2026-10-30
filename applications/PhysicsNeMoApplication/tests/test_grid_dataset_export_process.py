from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.processes.export import grid_dataset_export_process
from test_grid_bridge import CreateStructuredTetModelPart


class TestGridDatasetExportProcess(KratosUnittest.TestCase):
    def setUp(self):
        self.output_path = Path("test_grid_dataset_export")
        self.model = Kratos.Model()
        self.model_part = CreateStructuredTetModelPart(self.model, "Main", divisions=2)

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.output_path))

    def _CreateProcess(self, grid_shape=(4, 4, 4)):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "list_of_fields"  : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
                "grid_shape"      : [%d, %d, %d],
                "output_path"     : "%s",
                "output_interval" : 2
            }
        }""" % (grid_shape + (self.output_path,)))
        return grid_dataset_export_process.Factory(settings, self.model)

    def _SetField(self, value):
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, value)

    def test_ExportsGridSeriesOnInterval(self):
        process = self._CreateProcess()
        process.ExecuteInitialize()

        for step in range(1, 5):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            self.model_part.ProcessInfo[Kratos.TIME] = 0.5 * step
            self._SetField(float(step))
            process.ExecuteFinalizeSolutionStep()

        files = sorted(self.output_path.glob("*.npz"))
        self.assertEqual([f.name for f in files], ["grid_2.npz", "grid_4.npz"])

        with numpy.load(files[0]) as data:
            self.assertEqual(data["grid"].shape, (1, 4, 4, 4))
            self.assertEqual(data["grid"].dtype, numpy.float32)
            # constant field: every in-mesh lattice point carries it exactly
            self.assertTrue(numpy.allclose(data["grid"], 2.0))
            self.assertEqual(int(data["STEP"]), 2)
            self.assertAlmostEqual(float(data["TIME"]), 1.0)
            first_box = numpy.array(data["bounding_box"])
        with numpy.load(files[1]) as data:
            self.assertTrue(numpy.allclose(data["grid"], 4.0))
            # the bounding box is frozen after the first export
            numpy.testing.assert_allclose(numpy.array(data["bounding_box"]), first_box)

    def test_InvalidSettingsRaise(self):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "list_of_fields"  : [ { "variable_name" : "PRESSURE" } ],
                "grid_shape"      : [4, 4]
            }
        }""")
        with self.assertRaisesRegex(ValueError, "grid_shape"):
            grid_dataset_export_process.Factory(settings, self.model)

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "list_of_fields"  : [ { "variable_name" : "PRESSURE" } ],
                "output_interval" : 0
            }
        }""")
        with self.assertRaisesRegex(ValueError, "output_interval"):
            grid_dataset_export_process.Factory(settings, self.model)


if __name__ == '__main__':
    KratosUnittest.main()
