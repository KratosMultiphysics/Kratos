# --- Core Imports ---
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting

# --- HDF5 Imports ---
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application import core
from KratosMultiphysics.HDF5Application.core.utils import ParametersWrapper
from KratosMultiphysics.HDF5Application.core import controllers
from KratosMultiphysics.HDF5Application.core import operations
from KratosMultiphysics.HDF5Application.core import file_io

# --- STD Imports ---
import os
import pathlib
import typing
import functools
from unittest.mock import call, patch, MagicMock


test_file_path = pathlib.Path("kratos.h5")


def ensure_no_file(file_path: pathlib.Path) -> typing.Callable:
    """@brief Construct a decorator that deletes the specified file before and after the function is called."""
    def decorator(function: typing.Callable) -> typing.Callable:
        """@brief Delete a file before and after the function is called."""
        @functools.wraps(function)
        def wrapper(*args, **kwargs):
            if file_path.is_file():
                DeleteFileIfExisting(str(file_path))
            elif file_path.is_dir():
                raise FileExistsError(f"{file_path} is a directory")
            try:
                output = function(*args, **kwargs)
            finally:
                if file_path.is_file():
                    DeleteFileIfExisting(str(file_path))
                elif file_path.is_dir():
                    raise FileExistsError(f"{file_path} is a directory")
            return output
        return wrapper
    return decorator


def _SurrogateModelPart():
    model = KratosMultiphysics.Model()
    model_part = model.CreateModelPart("model_part")
    model_part.ProcessInfo[KratosMultiphysics.TIME] = 1.23456789
    model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = 0.1
    return model, model_part


class TestFileIO(KratosUnittest.TestCase):

    @staticmethod
    def _BuildTestFileIOObject(obj):
        class FilenameGetter(object):
            def Get(self, identifier):
                return identifier
        obj.filename_getter = FilenameGetter()
        obj.file_access_mode = 'exclusive'
        obj.file_driver = 'core'
        obj.echo_level = 0

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

    def test_FileIO_FileSettings(self):
        io = file_io._FileIO()
        self._BuildTestFileIOObject(io)
        settings = io._FileSettings('kratos.h5')
        self.assertEqual(settings['file_name'], 'kratos.h5')
        self.assertEqual(settings['file_access_mode'], 'exclusive')
        self.assertEqual(settings['file_driver'], 'core')
        self.assertEqual(settings['echo_level'], 0)

    def test_HDF5SerialFileIO_Creation(self):
        io = file_io._HDF5SerialFileIO()
        self._BuildTestFileIOObject(io)
        obj = io.Get('kratos.h5')
        self.assertIsInstance(obj, KratosHDF5.HDF5FileSerial)

    def test_SetDefaults(self):
        settings = ParametersWrapper()
        file_io._SetDefaults(settings)
        self.assertTrue(settings.Has('io_type'))
        self.assertTrue(settings.Has('file_name'))
        self.assertFalse(settings.Has('time_format'))
        self.assertTrue(settings.Has('file_access_mode'))
        self.assertTrue(settings.Has('file_driver'))
        self.assertTrue(settings.Has('echo_level'))
        self.assertEqual(settings['io_type'], 'serial_hdf5_file_io')
        self.assertEqual(settings['file_name'], 'kratos')
        self.assertEqual(settings['file_access_mode'], 'exclusive')
        if os.name == 'posix':
            self.assertEqual(settings['file_driver'], 'sec2')
        elif os.name == 'nt':
            self.assertEqual(settings['file_driver'], 'windows')
        self.assertEqual(settings['echo_level'], 0)

    def test_SetDefaults_NonTerminalTime(self):
        settings = ParametersWrapper()
        settings['file_name'] = 'kratos-<time>.h5'
        file_io._SetDefaults(settings)
        self.assertTrue(settings.Has('time_format'))
        self.assertEqual(settings['time_format'], '0.4f')

    def test_SetDefaults_ParallelIO(self):
        settings = ParametersWrapper()
        settings.AddEmptyValue('io_type').SetString('parallel_hdf5_file_io')
        file_io._SetDefaults(settings)
        self.assertEqual(settings['file_driver'], 'mpio')

    def test_GetIO_SerialIO(self):
        io = file_io._GetIO('serial_hdf5_file_io')
        self.assertIsInstance(io, file_io._HDF5SerialFileIO)

    def test_GetIO_ParallelIO(self):
        io = file_io._GetIO('parallel_hdf5_file_io')
        self.assertIsInstance(io, file_io._HDF5ParallelFileIO)

    def test_GetIO_MockIO(self):
        io = file_io._GetIO('mock_hdf5_file_io')
        self.assertIsInstance(io, file_io._HDF5MockFileIO)

    def test_GetIO_GarbageInput(self):
        with self.assertRaisesRegex(ValueError, r'"io_type" has invalid value "abcdefg"'):
            file_io._GetIO('abcdefg')

    @patch("KratosMultiphysics.FileNameDataCollector", autospec=True)
    def test_FilenameGetter_WithFileExtension(self, mock_class):
        settings = self._FilenameGetterSettings(file_name='kratos.h5')
        obj = file_io._FilenameGetter(settings)
        model, model_part = _SurrogateModelPart()
        obj.Get(model_part)
        mock_class.assert_called_once_with(model_part, 'kratos.h5', {'<time>':''})

    @patch("KratosMultiphysics.FileNameDataCollector", autospec=True)
    def test_FilenameGetter_WithoutFileExtension(self, mock_class):
        settings = self._FilenameGetterSettings(file_name='kratos')
        obj = file_io._FilenameGetter(settings)
        model, model_part = _SurrogateModelPart()
        obj.Get(model_part)
        mock_class.assert_called_once_with(model_part, 'kratos.h5', {'<time>':''})

    @patch("KratosMultiphysics.FileNameDataCollector", autospec=True)
    def test_FilenameGetter_TimeFormat(self, mock_class):
        settings = self._FilenameGetterSettings(time_format='0.4f')
        obj = file_io._FilenameGetter(settings)
        model, model_part = _SurrogateModelPart()
        obj.Get(model_part)
        mock_class.assert_called_once_with(model_part, 'kratos.h5', {'<time>':'0.4f'})

    @patch("KratosMultiphysics.FileNameDataCollector", autospec=True)
    def test_FilenameGetterWithDirectoryInitialization_WithoutDirectory(self, mock_class):
        mock_instance = mock_class.return_value
        mock_instance.GetFileName.return_value = 'kratos.h5'
        settings = self._FilenameGetterSettings()
        with patch('os.makedirs', autospec=True) as p:
            data_comm = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
            obj = file_io._FilenameGetterWithDirectoryInitialization(settings, data_comm)
            obj.Get(_SurrogateModelPart())
            self.assertEqual(p.call_count, 0)

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

    @patch("KratosMultiphysics.FileNameDataCollector", autospec=True)
    def test_FileIOMaxFilesToKeepExclusiveNoDeletion(self, mock_class):
        mock_class_instance = mock_class.return_value
        mock_class_instance.GetSortedFileNamesList.return_value = ['file_1', 'file_2', 'file_3', 'file_4', 'file_5']
        settings = self._FilenameGetterSettings(file_name='/foo/kratos.h5')
        settings['file_access_mode'] = 'exclusive'
        settings['max_files_to_keep'] = 4
        obj = file_io._FilenameGetter(settings)
        with patch("os.path.isdir", autospec=True) as mock_dir:
            mock_dir.return_value = True
            with patch("KratosMultiphysics.kratos_utilities.DeleteFileIfExisting", autospec=True) as p:
                obj.Get(_SurrogateModelPart())
                self.assertEqual(p.call_count, 0)
            self.assertEqual(mock_dir.call_args_list, [])

    @patch("KratosMultiphysics.FileNameDataCollector", autospec=True)
    def test_FileIOMaxFilesToKeepTruncateNoDeletion(self, mock_class):
        mock_class_instance = mock_class.return_value
        mock_class_instance.GetSortedFileNamesList.return_value = ['file_1', 'file_2', 'file_3']
        settings = self._FilenameGetterSettings(file_name='/foo/kratos.h5')
        settings['file_access_mode'] = 'truncate'
        settings['max_files_to_keep'] = 4
        obj = file_io._FilenameGetter(settings)
        with patch("pathlib.Path.parents", autospec=True) as mock_path:
            mock_is_dir = mock_path.__getitem__().is_dir
            mock_is_dir.return_value = True
            with patch("KratosMultiphysics.kratos_utilities.DeleteFileIfExisting", autospec=True) as p:
                obj.Get(_SurrogateModelPart())
                self.assertEqual(p.call_count, 0)
            self.assertEqual(mock_is_dir.call_count, 1)

    @patch("KratosMultiphysics.FileNameDataCollector", autospec=True)
    def test_FileIOMaxFilesToKeepTruncateDeletion(self, mock_class):
        mock_class_instance = mock_class.return_value
        mock_class_instance.GetSortedFileNamesList.return_value = ['file_1', 'file_2', 'file_3', 'file_4', 'file_5']
        settings = self._FilenameGetterSettings(file_name='/foo/kratos.h5')
        settings['file_access_mode'] = 'truncate'
        settings['max_files_to_keep'] = 4
        obj = file_io._FilenameGetter(settings)
        with patch("pathlib.Path.parents", autospec=True) as mock_path:
            mock_is_dir = mock_path.__getitem__().is_dir
            mock_is_dir.return_value = True
            with patch("KratosMultiphysics.kratos_utilities.DeleteFileIfExisting", autospec=True) as p:
                obj.Get(_SurrogateModelPart())
                self.assertEqual(p.call_args_list, [call('file_2'), call('file_3')])
            self.assertEqual(mock_is_dir.call_count, 1)

    def test_Create_Settings(self):
        settings = ParametersWrapper()
        data_comm = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
        file_io.Create(settings, data_comm)
        self.assertTrue(settings.Has('io_type'))
        self.assertTrue(settings.Has('file_name'))
        self.assertTrue(settings.Has('file_access_mode'))
        self.assertTrue(settings.Has('file_driver'))
        self.assertTrue(settings.Has('echo_level'))
        self.assertEqual(settings['io_type'], 'serial_hdf5_file_io')
        self.assertEqual(settings['file_access_mode'], 'exclusive')
        if os.name == 'posix':
            self.assertEqual(settings['file_driver'], 'sec2')
        elif os.name == 'nt':
            self.assertEqual(settings['file_driver'], 'windows')
        self.assertEqual(settings['echo_level'], 0)

    def test_Create_Attributes(self):
        settings = ParametersWrapper()
        data_comm = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
        io = file_io.Create(settings, data_comm)
        self.assertIsInstance(io, file_io._HDF5SerialFileIO)
        self.assertTrue(hasattr(io, 'filename_getter'))
        self.assertEqual(io.file_access_mode, 'exclusive')
        if os.name == 'posix':
            self.assertEqual(io.file_driver, 'sec2')
        elif os.name == 'nt':
            self.assertEqual(io.file_driver, 'windows')
        self.assertEqual(io.echo_level, 0)


class TestOperations(KratosUnittest.TestCase):

    @ensure_no_file(test_file_path)
    def test_CreateNonExistingOperation(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "abcdefg"
            }]
        }""")
        with self.assertRaises(ValueError):
            operations.Create(model, settings)

    @ensure_no_file(test_file_path)
    def test_Prefix_Literal(self):
        _, model_part = _SurrogateModelPart()
        prefix = operations.model_part.Prefix('/ModelData', model_part)
        self.assertEqual(prefix, '/ModelData')

    @ensure_no_file(test_file_path)
    def test_Prefix_NonTerminalTime(self):
        _, model_part = _SurrogateModelPart()
        prefix = operations.model_part.Prefix('/ModelData-<time>', model_part)
        self.assertEqual(prefix, '/ModelData-1.23456789')

    @ensure_no_file(test_file_path)
    def test_Prefix_FormattedNonTerminalTime(self):
        _, model_part = _SurrogateModelPart()
        prefix = operations.model_part.Prefix(
            '/ModelData-<time>', model_part, '0.2f')
        self.assertEqual(prefix, '/ModelData-1.23')

    @ensure_no_file(test_file_path)
    def test_Prefix_NonTerminalIdentifier(self):
        _, model_part = _SurrogateModelPart()
        prefix = operations.model_part.Prefix(
            '/<model_part_name>-<time>', model_part)
        self.assertEqual(prefix, '/model_part-1.23456789')

    @ensure_no_file(test_file_path)
    def test_ModelPartOutput(self):
        model, model_part = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "model_part_output"
            }]
        }""")
        operation_settings = settings["list_of_operations"][0]
        model_part_output = operations.Create(model, settings)
        self.assertTrue(operation_settings.Has('prefix'))
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ModelPartIO', autospec=True) as p:
            model_part_io = p.return_value
            model_part_output.Execute()
            model_part_io.WriteModelPart.assert_called_once_with(model_part)

    @ensure_no_file(test_file_path)
    def test_ModelPartOutput_NonTerminalPrefix(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "model_part_output",
                "prefix": "/ModelData/<model_part_name>/<time>",
                "time_format": "0.2f"
            }]
        }""")
        model_part_output = operations.Create(model, settings)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ModelPartIO', autospec=True) as p:
            model_part_output.Execute()
            args, _ = p.call_args
            self.assertEqual(args[1], '/ModelData/model_part/1.23')

    @ensure_no_file(test_file_path)
    def test_NodalSolutionStepDataOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "nodal_solution_step_data_output"
            }]
        }""")
        nodal_solution_step_data_output = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepDataIO', autospec=True) as p:
            nodal_solution_step_data_io = p.return_value
            nodal_solution_step_data_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_solution_step_data_io.WriteNodalResults.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_NodalSolutionStepDataInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "nodal_solution_step_data_input"
            }]
        }""")
        nodal_solution_step_data_input = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepDataIO', autospec=True) as p:
            nodal_solution_step_data_io = p.return_value
            nodal_solution_step_data_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_solution_step_data_io.ReadNodalResults.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_NodalDataValueOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "nodal_data_value_output"
            }]
        }""")
        nodal_data_value_output = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalDataValueIO', autospec=True) as p:
            nodal_data_value_io = p.return_value
            nodal_data_value_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_data_value_io.WriteNodalResults.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_NodalFlagValueOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "nodal_flag_value_output"
            }]
        }""")
        nodal_flag_value_output = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalFlagValueIO', autospec=True) as p:
            nodal_flag_value_io = p.return_value
            nodal_flag_value_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_flag_value_io.WriteNodalFlags.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_NodalDataValueInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "nodal_data_value_input"
            }]
        }""")
        nodal_data_value_input = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalDataValueIO', autospec=True) as p:
            nodal_data_value_io = p.return_value
            nodal_data_value_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_data_value_io.ReadNodalResults.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_NodalFlagValueInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "nodal_flag_value_input"
            }]
        }""")
        nodal_flag_value_input = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalFlagValueIO', autospec=True) as p:
            nodal_flag_value_io = p.return_value
            nodal_flag_value_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_flag_value_io.ReadNodalFlags.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ElementDataValueOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "element_data_value_output"
            }]
        }""")
        element_data_value_output = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementDataValueIO', autospec=True) as p:
            element_data_value_io = p.return_value
            element_data_value_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(element_data_value_io.WriteElementResults.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ElementFlagValueOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "element_flag_value_output"
            }]
        }""")
        element_flag_value_output = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementFlagValueIO', autospec=True) as p:
            element_flag_value_io = p.return_value
            element_flag_value_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(element_flag_value_io.WriteElementFlags.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ElementDataValueInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "element_data_value_input"
            }]
        }""")
        element_data_value_input = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementDataValueIO', autospec=True) as p:
            element_data_value_io = p.return_value
            element_data_value_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(element_data_value_io.ReadElementResults.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ElementFlagValueInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "element_flag_value_input"
            }]
        }""")
        element_flag_value_input = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementFlagValueIO', autospec=True) as p:
            element_flag_value_io = p.return_value
            element_flag_value_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(element_flag_value_io.ReadElementFlags.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ConditionDataValueOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "condition_data_value_output"
            }]
        }""")
        condition_data_value_output = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionDataValueIO', autospec=True) as p:
            condition_data_value_io = p.return_value
            condition_data_value_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(condition_data_value_io.WriteConditionResults.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ConditionFlagValueOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "condition_flag_value_output"
            }]
        }""")
        condition_flag_value_output = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionFlagValueIO', autospec=True) as p:
            condition_flag_value_io = p.return_value
            condition_flag_value_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(condition_flag_value_io.WriteConditionFlags.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ConditionDataValueInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "condition_data_value_input"
            }]
        }""")
        condition_data_value_input = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionDataValueIO', autospec=True) as p:
            condition_data_value_io = p.return_value
            condition_data_value_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(condition_data_value_io.ReadConditionResults.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ConditionFlagValueInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "condition_flag_value_input"
            }]
        }""")
        condition_flag_value_input = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionFlagValueIO', autospec=True) as p:
            condition_flag_value_io = p.return_value
            condition_flag_value_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(condition_flag_value_io.ReadConditionFlags.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_PrimalBossakOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "primal_bossak_output"
            }]
        }""")
        primal_bossak_output = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepBossakIO', autospec=True) as p:
            nodal_solution_step_bossak_io = p.return_value
            primal_bossak_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_solution_step_bossak_io.WriteNodalResults.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_PrimalBossakInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "primal_bossak_input"
            }]
        }""")
        primal_bossak_input = operations.Create(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepBossakIO', autospec=True) as p:
            nodal_solution_step_bossak_io = p.return_value
            primal_bossak_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_solution_step_bossak_io.ReadNodalResults.call_count, 1)


class TestControllers(KratosUnittest.TestCase):

    def test_CreateNonExistingController(self):
        settings = ParametersWrapper()
        settings['controller_type'] = 'abcdefg'
        with self.assertRaisesRegex(ValueError, r'"controller_type" has invalid value "abcdefg"'):
            controllers.Factory(MagicMock(), MagicMock(), settings.Get())

    @patch('KratosMultiphysics.FileNameDataCollector', autospec=True)
    def test_DefaultController(self, mock_class):
        mock_instance = mock_class.return_value
        mock_instance.GetFileName.return_value = 'kratos.h5'
        _, model_part = _SurrogateModelPart()
        controller_settings = ParametersWrapper()
        with patch('KratosMultiphysics.HDF5Application.core.file_io.KratosHDF5.HDF5FileSerial', autospec=True):
            operation = MagicMock(spec = operations.AggregateOperation)
            controller = controllers.Factory(
                model_part, operation, controller_settings.Get())
            self.assertTrue(controller_settings.Has('controller_type'))
            self.assertEqual(controller_settings['controller_type'], 'default_controller')
            for _ in range(10):
                controller()
            self.assertEqual(operation.Execute.call_count, 10)

    def test_TemporalController_CreateWithDefaults(self):
        _, model_part = _SurrogateModelPart()
        operation = MagicMock(spec = operations.AggregateOperation)
        controller_settings = ParametersWrapper()
        controller_settings['controller_type'] = 'temporal_controller'
        controller = controllers.Factory(model_part, operation, controller_settings.Get())
        for _ in range(1, 11):
            model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
            model_part.ProcessInfo[KratosMultiphysics.TIME] += model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
            controller()
        self.assertEqual(operation.Execute.call_count, 10)

    def test_TemporalController_CreateWithParameters(self):
        _, model_part = _SurrogateModelPart()
        operation = MagicMock(spec = operations.AggregateOperation)
        controller_settings = ParametersWrapper()
        controller_settings['controller_type'] = 'temporal_controller'
        controller_settings['time_frequency'] = 2.0
        controller_settings['step_frequency'] = 3
        controller = controllers.Factory(model_part, operation, controller_settings.Get())
        for _ in range(1, 31):
            model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
            model_part.ProcessInfo[KratosMultiphysics.TIME] += model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
            controller()
        self.assertEqual(operation.Execute.call_count, 10)

    def test_TemporalController_StepFrequency(self):
        _, model_part = _SurrogateModelPart()
        controller_settings = ParametersWrapper()
        controller_settings['step_frequency'] = 2
        controller_settings['controller_type'] = 'temporal_controller'
        with patch('KratosMultiphysics.HDF5Application.core.file_io._HDF5SerialFileIO', autospec=True):
            operation = MagicMock(spec = operations.AggregateOperation)
            controller = controllers.Factory(
                model_part, operation, controller_settings.Get())
            for _ in range(10):
                model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
                model_part.ProcessInfo[KratosMultiphysics.TIME] += model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
                controller()
            self.assertEqual(operation.Execute.call_count, 5)

    def test_TemporalController_TimeFrequency(self):
        _, model_part = _SurrogateModelPart()
        controller_settings = ParametersWrapper()
        controller_settings['step_frequency'] = 100
        controller_settings['time_frequency'] = 0.5
        controller_settings['controller_type'] = 'temporal_controller'
        with patch('KratosMultiphysics.HDF5Application.core.file_io._HDF5SerialFileIO', autospec=True):
            operation = MagicMock(spec = operations.AggregateOperation)
            controller = controllers.Factory(
                model_part, operation, controller_settings.Get())
            for _ in range(10):
                model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
                model_part.ProcessInfo[KratosMultiphysics.TIME] += model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
                controller()
            self.assertEqual(operation.Execute.call_count, 2)

    def test_TemporalController_NearlyTheSameTimeFrequency(self):
        _, model_part = _SurrogateModelPart()
        with patch('KratosMultiphysics.HDF5Application.core.file_io._HDF5SerialFileIO', autospec=True):
            controller_settings = ParametersWrapper()
            controller_settings['step_frequency'] = 100
            controller_settings['time_frequency'] = 0.2000001
            controller_settings['controller_type'] = 'temporal_controller'
            operation = MagicMock(spec = operations.AggregateOperation)
            controller = controllers.Factory(model_part, operation, controller_settings.Get())
            for _ in range(10):
                model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
                model_part.ProcessInfo[KratosMultiphysics.TIME] += model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
                controller()
            self.assertEqual(operation.Execute.call_count, 5)

    @patch('KratosMultiphysics.FileNameDataCollector', autospec=True)
    def test_TemporalController_OperationCall(self, mock_class):
        mock_instance = mock_class.return_value
        mock_instance.GetFileName.return_value = 'kratos.h5'
        _, model_part = _SurrogateModelPart()
        controller_settings = KratosMultiphysics.Parameters("""{
            "controller_type" : "temporal_controller"
        }""")
        operation = MagicMock(spec = operations.AggregateOperation)
        controller = controllers.Factory(
            model_part, operation, controller_settings)
        with patch('KratosMultiphysics.HDF5Application.core.file_io.KratosHDF5.HDF5FileSerial', autospec=True):
            for _ in range(10):
                model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
                model_part.ProcessInfo[KratosMultiphysics.TIME] += model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
                controller()
            self.assertEqual(operation.Execute.call_count, 10)


class TestFactory(KratosUnittest.TestCase):

    def test_NonArraySettings(self):
        model = KratosMultiphysics.Model()
        settings = KratosMultiphysics.Parameters()
        with self.assertRaisesRegex(ValueError, r'Expected settings as an array'):
            core.Factory(settings, model, KratosMultiphysics.Process)

    def test_EmptyArraySettings(self):
        model = KratosMultiphysics.Model()
        settings = KratosMultiphysics.Parameters('''
            {
                "list_of_controllers" : []
            }
            ''')
        settings = ParametersWrapper(settings)
        with self.assertRaises(Exception):
            core.Factory(settings['list_of_controllers'], model, KratosMultiphysics.Process)

    def test_DefaultSettings(self):
        model = KratosMultiphysics.Model()
        model.CreateModelPart('test')
        parent_settings = KratosMultiphysics.Parameters('''{
                "list_of_controllers" : [{
                        "model_part_name" : "test"
                }]
            }''')
        parent_settings = ParametersWrapper(parent_settings)
        core.Factory(parent_settings['list_of_controllers'], model, KratosMultiphysics.Process)
        settings = parent_settings['list_of_controllers'][0]
        self.assertTrue(settings.Has('model_part_name'))
        self.assertTrue(settings.Has('process_step'))
        self.assertTrue(settings.Has('controller_settings'))
        self.assertTrue(settings.Has('io_settings'))
        self.assertTrue(settings.Has('list_of_operations'))
        self.assertTrue(settings['list_of_operations'].IsArray())
        self.assertEqual(settings['list_of_operations'].size(), 0)

    def test_DefaultProcess(self):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart('test')
        parent_settings = KratosMultiphysics.Parameters('''{
                "list_of_controllers" : [{
                    "model_part_name" : "test",
                    "list_of_operations" : [{
                        "operation_type" : "model_part_output"
                    }]
                }]
            }''')
        parent_settings = ParametersWrapper(parent_settings)
        process = core.Factory(parent_settings['list_of_controllers'], model, KratosMultiphysics.Process)
        with patch('KratosMultiphysics.HDF5Application.core.file_io.KratosHDF5.HDF5FileSerial', autospec=True) as MockedFileSerial:
            with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ModelPartIO', autospec=True) as MockedModelPartIO:
                process.ExecuteInitialize()
                model_part_io = MockedModelPartIO.return_value
                model_part_io.WriteModelPart.assert_called_once_with(model_part)


    def test_TemporalController_PrintOutput(self):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("test")
        parameters = ParametersWrapper("""[{
            "model_part_name" : "test",
            "process_step" : "finalize_solution_step",
            "controller_settings" : {
                "controller_type" : "temporal_controller",
                "time_frequency" : 5.0,
                "step_frequency" : 5
            },
            "list_of_operations" : [{"operation_type" : "model_part_output"}]
        }]""")
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ModelPartIO', autospec=True) as MockedModelPartIO:
            with patch("KratosMultiphysics.HDF5Application.core.file_io._HDF5SerialFileIO", autospec = True) as MockedSerialFileIO:
                process = core.Factory(parameters, model, KratosMultiphysics.OutputProcess)

                process.ExecuteInitialize()
                process.ExecuteBeforeSolutionLoop()

                MockedModelPartIO.reset_mock() # discard initial calls

                # No writes should be performed here
                for step in range(4):
                    model_part.CloneTimeStep(step)
                    model_part.ProcessInfo[KratosMultiphysics.STEP] = step
                    process.ExecuteFinalizeSolutionStep()
                process.ExecuteFinalize()
                MockedModelPartIO.assert_not_called()

                # This should trigger a single write
                process.PrintOutput()
                MockedModelPartIO.assert_called_once()


    def test_TemporalController_OutputStep(self):
        model = KratosMultiphysics.Model()
        model.CreateModelPart("test")
        parameters = ParametersWrapper("""[{
            "model_part_name" : "test",
            "process_step" : "output",
            "controller_settings" : {
                "controller_type" : "temporal_controller"
            },
            "list_of_operations" : []
        }]""")
        with patch("KratosMultiphysics.HDF5Application.core.file_io._HDF5SerialFileIO", autospec = True):
            with self.assertRaises(TypeError):
                core.Factory(parameters, model, KratosMultiphysics.OutputProcess)


class TestParametersWrapper(KratosUnittest.TestCase):

    def setUp(self):
        self.set_params = KratosMultiphysics.Parameters()
        self.get_params = KratosMultiphysics.Parameters(
            '''
            {
            "string_value" : "abc",
            "int_value": 1,
            "double_value": 1.5,
            "bool_value": true,
            "parameters" : {
                "double_value": 3.1
            },
            "array_of_double_values": [1.1, 2.2, 3.3],
            "array_of_parameters": [{
                "int_value": 2
            },{
                "double_value": 2.7
            }]
            }
            '''
        )

    def test_get_string(self):
        settings = core.ParametersWrapper(self.get_params)
        self.assertEqual(settings['string_value'], 'abc')

    def test_get_int(self):
        settings = core.ParametersWrapper(self.get_params)
        self.assertEqual(settings['int_value'], 1)

    def test_get_double(self):
        settings = core.ParametersWrapper(self.get_params)
        self.assertAlmostEqual(settings['double_value'], 1.5)

    def test_get_bool(self):
        settings = core.ParametersWrapper(self.get_params)
        self.assertEqual(settings['bool_value'], True)

    def test_get_parameters(self):
        settings = core.ParametersWrapper(self.get_params)
        self.assertAlmostEqual(settings['parameters']['double_value'], 3.1)

    def test_get_array_of_values(self):
        settings = core.ParametersWrapper(self.get_params)
        self.assertAlmostEqual(settings['array_of_double_values'][0], 1.1)
        self.assertAlmostEqual(settings['array_of_double_values'][1], 2.2)
        self.assertAlmostEqual(settings['array_of_double_values'][2], 3.3)

    def test_get_array_of_parameters(self):
        settings = core.ParametersWrapper(self.get_params)
        self.assertEqual(settings['array_of_parameters'][0]['int_value'], 2)
        self.assertEqual(settings['array_of_parameters'][1]['double_value'], 2.7)

    def test_set_string(self):
        settings = core.ParametersWrapper(self.set_params)
        settings['string_value'] = 'abc'
        self.assertEqual(settings['string_value'], 'abc')

    def test_set_int(self):
        settings = core.ParametersWrapper(self.set_params)
        settings['int_value'] = 1
        self.assertEqual(settings['int_value'], 1)

    def test_set_double(self):
        settings = core.ParametersWrapper(self.set_params)
        settings['double_value'] = 1.5
        self.assertEqual(settings['double_value'], 1.5)

    def test_set_bool(self):
        settings = core.ParametersWrapper(self.set_params)
        settings['bool_value'] = True
        self.assertEqual(settings['bool_value'], True)

    def test_set_parameters(self):
        settings = core.ParametersWrapper(self.set_params)
        settings['parameters'] = KratosMultiphysics.Parameters(
            '''
            {
            "bool_value": false
            }
            '''
        )
        self.assertAlmostEqual(settings['parameters']['bool_value'], False)

    def test_set_array_of_values(self):
        settings = core.ParametersWrapper(self.set_params)
        settings['array_of_bool_values'] = [True, False]
        self.assertAlmostEqual(settings['array_of_bool_values'][0], True)
        self.assertAlmostEqual(settings['array_of_bool_values'][1], False)

    def test_set_array_of_parameters(self):
        settings = core.ParametersWrapper(self.set_params)
        settings['array_of_parameters'] = [
            self.get_params['array_of_parameters'][0],
            self.get_params['array_of_parameters'][1]
        ]
        self.assertEqual(settings['array_of_parameters'][0]['int_value'], 2)
        self.assertEqual(settings['array_of_parameters'][1]['double_value'], 2.7)

    def test_array_keys(self):
        settings = ParametersWrapper(self.get_params)
        count = 0
        for k in settings['array_of_double_values']:
            self.assertEqual(k, count)
            count += 1

    def test_nonarray_keys(self):
        settings = ParametersWrapper(self.get_params)
        for k in settings['parameters']:
            self.assertEqual(k, 'double_value')


if __name__ == "__main__":
    KratosUnittest.main()
