import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application import core
from KratosMultiphysics.HDF5Application.core import operations, file_io
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.HDF5Application.core.utils import ParametersWrapper
from unittest.mock import patch, MagicMock
import test_hdf5_core

def _SurrogateModelPart():
    model_part = MagicMock(spec=KratosMultiphysics.ModelPart)
    model_part.ProcessInfo = {KratosMultiphysics.TIME: 1.23456789,
                              KratosMultiphysics.DELTA_TIME: 0.1}
    model_part.Name = 'model_part'
    return model_part

class TestFileIO(KratosUnittest.TestCase):

    def test_HDF5ParallelFileIO_Creation(self):
        io = file_io._HDF5ParallelFileIO()
        test_hdf5_core.TestFileIO._BuildTestFileIOObject(io)
        obj = io.Get('kratos.h5')
        self.assertIsInstance(obj, KratosHDF5.HDF5FileParallel)

    @patch("KratosMultiphysics.FileNameDataCollector", autospec=True)
    def test_FilenameGetterWithDirectoryInitialization_DirectoryExists(self, mock_class):
        mock_instance = mock_class.return_value
        mock_instance.GetFileName.return_value = '/foo/kratos.h5'
        settings = self._FilenameGetterSettings(file_name='/foo/kratos.h5')
        patcher = patch('KratosMultiphysics.FilesystemExtensions.MPISafeCreateDirectories', autospec=True)
        makedirs = patcher.start()
        data_comm = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
        obj = file_io._FilenameGetterWithDirectoryInitialization(settings, data_comm)
        obj.Get(_SurrogateModelPart())
        makedirs.assert_called_once_with('/foo')
        patcher.stop()

    @patch("KratosMultiphysics.FileNameDataCollector", autospec=True)
    def test_FilenameGetterWithDirectoryInitialization_DirectoryDoesNotExist(self, mock_class):
        mock_instance = mock_class.return_value
        mock_instance.GetFileName.return_value = '/foo/kratos.h5'
        settings = self._FilenameGetterSettings(file_name='/foo/kratos.h5')
        patcher = patch('KratosMultiphysics.FilesystemExtensions.MPISafeCreateDirectories', autospec=True)
        makedirs = patcher.start()
        data_comm = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
        obj = file_io._FilenameGetterWithDirectoryInitialization(settings, data_comm)
        obj.Get(_SurrogateModelPart())
        makedirs.assert_called_once_with('/foo')
        patcher.stop()

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
        settings = ParametersWrapper()
        settings['operation_type'] = 'partitioned_model_part_output'
        partitioned_model_part_output = operations.Create(settings)
        self.assertTrue(settings.Has('operation_type'))
        self.assertTrue(settings.Has('prefix'))
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5PartitionedModelPartIO') as p:
            partitioned_model_part_io = p.return_value
            model_part = test_hdf5_core._SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileParallel)
            partitioned_model_part_output(model_part, hdf5_file)
            p.assert_called_once_with(hdf5_file, '/ModelData')
            partitioned_model_part_io.WriteModelPart.assert_called_once_with(
                model_part)

    def test_PartitionedModelPartOutput_NonTerminalPrefix(self):
        settings = ParametersWrapper('''
            {
                "operation_type": "partitioned_model_part_output",
                "prefix": "/ModelData/<model_part_name>/<time>",
                "time_format": "0.2f"
            }
            ''')
        partitioned_model_part_output = operations.Create(settings)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5PartitionedModelPartIO', autospec=True) as p:
            model_part = test_hdf5_core._SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileParallel)
            partitioned_model_part_output(model_part, hdf5_file)
            args, _ = p.call_args
            self.assertEqual(args[1], '/ModelData/model_part/1.23')

if __name__ == "__main__":
    KratosUnittest.main()
