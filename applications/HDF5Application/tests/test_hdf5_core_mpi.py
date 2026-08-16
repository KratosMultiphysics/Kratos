# --- Core Imports ---
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import ReadModelPart

# --- HDF5 Imports ---
import KratosMultiphysics.HDF5Application.core.operations as operations

# --- STD Imports ---
from unittest.mock import patch
import pathlib


def _SurrogateModelPart():
    # Get input model part path
    script_directory      = pathlib.Path(__file__).absolute().parent
    kratos_root_directory = script_directory.parent.parent.parent
    test_input_directory  = kratos_root_directory / "kratos" / "tests" / "auxiliar_files_for_python_unittest"
    test_file_stem        = test_input_directory / "mdpa_files" / "test_processes"
    test_file_path        = pathlib.Path(str(test_file_stem) + ".mdpa")

    if not test_file_path.is_file():
        raise FileNotFoundError(f"Test file not found: {test_file_path}")

    # Construct model and model part
    model = KratosMultiphysics.Model()
    model_part = model.CreateModelPart("model_part")

    # Read and set model part
    ReadModelPart(str(test_file_stem), model_part)
    model_part.ProcessInfo[KratosMultiphysics.TIME] = 1.23456789
    model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = 0.1
    return model, model_part


class TestOperations(KratosUnittest.TestCase):

    def test_PartitionedModelPartOutput(self):
        model, model_part = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "partitioned_model_part_output"
            }]
        }""")
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(model_part.IsDistributed())
        with patch("KratosMultiphysics.HDF5Application.HDF5File"):
            partitioned_model_part_output = operations.CreateAggregatedOperation(model, settings)
            self.assertTrue(operation_settings.Has('prefix'))
            with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5PartitionedModelPartIO') as p:
                partitioned_model_part_io = p.return_value
                partitioned_model_part_output.Execute()
                partitioned_model_part_io.WriteModelPart.assert_called_once_with(model_part)

    def test_PartitionedModelPartOutput_NonTerminalPrefix(self):
        model, model_part = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type": "partitioned_model_part_output",
                "prefix": "/ModelData/<model_part_name>/<time>",
                "time_format": "0.2f"
            }]
        }""")
        self.assertTrue(model_part.IsDistributed())
        with patch("KratosMultiphysics.HDF5Application.HDF5File"):
            partitioned_model_part_output = operations.CreateAggregatedOperation(model, settings)
            with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5PartitionedModelPartIO', autospec=True) as p:
                partitioned_model_part_output.Execute()
                args, _ = p.call_args
                self.assertEqual(args[0]["prefix"].GetString(), '/ModelData/model_part/1.23')

if __name__ == "__main__":
    KratosUnittest.main()
