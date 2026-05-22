"""
Serial tests for KaHIPPartitioningModeler.

Tests the modeler in single-process mode (no MPI), verifying that it correctly
reads a mesh file into a model part when perform_partitioning=false.
"""

import pathlib

import KratosMultiphysics as KM
from KratosMultiphysics import (
    KratosUnittest,
    Model,
    Logger,
    Parameters,
    PARTITION_INDEX,
)
from KratosMultiphysics.KaHIPApplication import KaHIPPartitioningModeler


WORK_DIR = pathlib.Path(__file__).parent.absolute()


class TestKaHIPModeler(KratosUnittest.TestCase):
    """Tests for KaHIPPartitioningModeler in serial (no-partition) mode."""

    # ── Basic construction ────────────────────────────────────────────────────

    def test_default_parameters(self):
        """GetDefaultParameters() must return a valid Parameters object."""
        model = Model()
        settings = Parameters("""{
            "model_part_name": "Test",
            "input_filename": "dummy"
        }""")
        modeler = KaHIPPartitioningModeler(model, settings)
        default_params = modeler.GetDefaultParameters()

        self.assertEqual(default_params["model_part_name"].GetString(), "Main")
        self.assertEqual(default_params["number_of_partitions"].GetInt(), 0)
        self.assertFalse(default_params["partition_in_memory"].GetBool())
        self.assertTrue(default_params["perform_partitioning"].GetBool())
        self.assertEqual(
            default_params["kahip_settings"]["preconfiguration"].GetString(), "eco")

    def test_info_string(self):
        """Info() must return the expected class name."""
        model = Model()
        settings = Parameters("""{
            "model_part_name": "Test",
            "input_filename": "dummy"
        }""")
        modeler = KaHIPPartitioningModeler(model, settings)
        self.assertEqual(modeler.Info(), "KaHIPPartitioningModeler")

    # ── Serial read (no partitioning) ─────────────────────────────────────────

    def test_no_partition_reads_model_part(self):
        """perform_partitioning=false should read the mesh directly without KaHIP."""
        file_name = str(WORK_DIR / "quads")
        model = Model()
        model_part = model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)

        settings = Parameters("""{
            "model_part_name": "Main",
            "input_filename": \"""" + file_name + """\",
            "perform_partitioning": false,
            "skip_timer": true,
            "ignore_variables_not_in_solution_step_data": true
        }""")

        modeler = KaHIPPartitioningModeler(model, settings)
        modeler.SetupModelPart()

        self.assertGreater(model_part.NumberOfNodes(), 0)
        self.assertGreater(model_part.NumberOfElements(), 0)

    def test_no_partition_cube(self):
        """Read the cube.mdpa mesh directly without partitioning."""
        file_name = str(WORK_DIR / "cube")
        model = Model()
        model_part = model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)

        settings = Parameters("""{
            "model_part_name": "Main",
            "input_filename": \"""" + file_name + """\",
            "perform_partitioning": false,
            "skip_timer": true,
            "ignore_variables_not_in_solution_step_data": true
        }""")

        modeler = KaHIPPartitioningModeler(model, settings)
        modeler.SetupModelPart()

        self.assertEqual(model_part.NumberOfNodes(), 413)
        self.assertEqual(model_part.NumberOfElements(), 1191)
        self.assertEqual(model_part.NumberOfConditions(), 780)

    def test_missing_filename_raises(self):
        """Empty input_filename must raise an error."""
        model = Model()
        settings = Parameters("""{
            "model_part_name": "Main",
            "input_filename": ""
        }""")
        modeler = KaHIPPartitioningModeler(model, settings)
        with self.assertRaises(Exception):
            modeler.SetupModelPart()

    # ── Single-process MPI path (mpi_size == 1 → no partitioning) ────────────

    def test_single_process_creates_model_part_if_missing(self):
        """Model part is created automatically when not present in the Model."""
        file_name = str(WORK_DIR / "quads")
        model = Model()  # model part "AutoCreated" does not exist yet

        settings = Parameters("""{
            "model_part_name": "AutoCreated",
            "input_filename": \"""" + file_name + """\",
            "perform_partitioning": false,
            "skip_timer": true,
            "ignore_variables_not_in_solution_step_data": true
        }""")

        modeler = KaHIPPartitioningModeler(model, settings)
        modeler.SetupModelPart()

        self.assertTrue(model.HasModelPart("AutoCreated"))
        mp = model.GetModelPart("AutoCreated")
        self.assertGreater(mp.NumberOfNodes(), 0)


if __name__ == '__main__':
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
