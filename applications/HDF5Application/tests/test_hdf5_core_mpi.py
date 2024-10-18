# --- Core Imports ---
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import ReadModelPart

# --- HDF5 Imports ---
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application.core import operations, file_io
from KratosMultiphysics.HDF5Application.core.utils import ParametersWrapper
import test_hdf5_core

# --- STD Imports ---
from unittest.mock import patch
import pathlib


def _SurrogateModelPart():
    # Get input model part path
    script_directory      = pathlib.Path(__file__).absolute().parent
    kratos_root_directory = script_directory.parent.parent.parent
    test_input_directory  = kratos_root_directory / "kratos" / "tests" / "test_files"
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


class TestFileIO(KratosUnittest.TestCase):

    def test_HDF5ParallelFileIO_Creation(self):
        io = file_io._HDF5ParallelFileIO()
        test_hdf5_core.TestFileIO._BuildTestFileIOObject(io)
        obj = io.Get('kratos.h5')
        self.assertIsInstance(obj, KratosHDF5.HDF5File)

    @patch("KratosMultiphysics.FileNameDataCollector", autospec=True)
    def test_FilenameGetterWithDirectoryInitialization_DirectoryExists(self, mock_class):
        mock_instance = mock_class.return_value
        mock_instance.GetFileName.return_value = '/foo/kratos.h5'
        settings = self._FilenameGetterSettings(file_name='/foo/kratos.h5')
        with patch('KratosMultiphysics.FilesystemExtensions.MPISafeCreateDirectories', autospec=True) as makedirs:
            data_comm = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
            obj = file_io._FilenameGetterWithDirectoryInitialization(settings, data_comm)
            obj.Get(_SurrogateModelPart()[1])
            makedirs.assert_called_once_with('/foo')

    @patch("KratosMultiphysics.FileNameDataCollector", autospec=True)
    def test_FilenameGetterWithDirectoryInitialization_DirectoryDoesNotExist(self, mock_class):
        mock_instance = mock_class.return_value
        mock_instance.GetFileName.return_value = '/foo/kratos.h5'
        settings = self._FilenameGetterSettings(file_name='/foo/kratos.h5')
        with patch('KratosMultiphysics.FilesystemExtensions.MPISafeCreateDirectories', autospec=True) as makedirs:
            data_comm = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
            obj = file_io._FilenameGetterWithDirectoryInitialization(settings, data_comm)
            obj.Get(_SurrogateModelPart()[1])
            makedirs.assert_called_once_with('/foo')

    @staticmethod
    def _FilenameGetterSettings(**kwargs):
        settings = ParametersWrapper()
        if 'file_name' in kwargs:
            settings['file_name'] = kwargs['file_name']
        else:
            settings['file_name'] = 'kratos.h5'
        if 'time_format' in kwargs:
            settings['time_format'] = kwargs['time_format']
        return settings


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
            partitioned_model_part_output = operations.Create(model, settings)
            self.assertTrue(operation_settings.Has('operation_type'))
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
            partitioned_model_part_output = operations.Create(model, settings)
            with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5PartitionedModelPartIO', autospec=True) as p:
                partitioned_model_part_output.Execute()
                args, _ = p.call_args
                self.assertEqual(args[1], '/ModelData/model_part/1.23')

if __name__ == "__main__":
    KratosUnittest.main()
