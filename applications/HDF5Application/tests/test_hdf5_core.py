import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application import core
from KratosMultiphysics.HDF5Application.core.utils import ParametersWrapper
from KratosMultiphysics.HDF5Application.core import controllers
from KratosMultiphysics.HDF5Application.core import operations
from KratosMultiphysics.HDF5Application.core import file_io
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
from unittest.mock import patch, MagicMock


def _SurrogateModelPart():
    model_part = MagicMock(spec=KratosMultiphysics.ModelPart)
    model_part.ProcessInfo = {KratosMultiphysics.TIME: 1.23456789,
                              KratosMultiphysics.DELTA_TIME: 0.1}
    model_part.Name = 'model_part'
    return model_part


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

    def test_FilenameGetter_WithFileExtension(self):
        settings = self._FilenameGetterSettings(file_name='kratos.h5')
        obj = file_io._FilenameGetter(settings)
        self.assertEqual(obj.Get(), 'kratos.h5')

    def test_FilenameGetter_WithoutFileExtension(self):
        settings = self._FilenameGetterSettings(file_name='kratos')
        obj = file_io._FilenameGetter(settings)
        self.assertEqual(obj.Get(), 'kratos.h5')

    def test_FilenameGetter_TimeFormat(self):
        settings = self._FilenameGetterSettings(time_format='0.4f')
        obj = file_io._FilenameGetter(settings)
        self.assertEqual(obj.Get(_SurrogateModelPart()), 'kratos.h5')

    def test_FilenameGetter_NonTerminalTime(self):
        settings = self._FilenameGetterSettings(file_name='kratos-<time>.h5')
        obj = file_io._FilenameGetter(settings)
        self.assertEqual(obj.Get(_SurrogateModelPart()),
                         'kratos-1.23456789.h5')

    def test_FilenameGetter_FormattedNonTerminalTime(self):
        settings = self._FilenameGetterSettings(
            file_name='kratos-<time>.h5', time_format='0.2f')
        obj = file_io._FilenameGetter(settings)
        self.assertEqual(obj.Get(_SurrogateModelPart()), 'kratos-1.23.h5')

    def test_FilenameGetter_NonTerminalIdentifier(self):
        settings = self._FilenameGetterSettings(
            file_name='<model_part_name>-<time>.h5', time_format='0.2f')
        obj = file_io._FilenameGetter(settings)
        self.assertEqual(obj.Get(_SurrogateModelPart()), 'model_part-1.23.h5')

    def test_FilenameGetterWithDirectoryInitialization_WithoutDirectory(self):
        settings = self._FilenameGetterSettings()
        with patch('os.makedirs', autospec=True) as p:
            obj = file_io._FilenameGetterWithDirectoryInitialization(settings)
            obj.Get(_SurrogateModelPart())
            self.assertEqual(p.call_count, 0)

    def test_FilenameGetterWithDirectoryInitialization_DirectoryExists(self):
        settings = self._FilenameGetterSettings(file_name='/foo/kratos.h5')
        patcher1 = patch('os.path.exists', autospec=True)
        patcher2 = patch('os.makedirs', autospec=True)
        pathexists = patcher1.start()
        makedirs = patcher2.start()
        pathexists.return_value = True
        obj = file_io._FilenameGetterWithDirectoryInitialization(settings)
        obj.Get(_SurrogateModelPart())
        pathexists.assert_called_once_with('/foo')
        self.assertEqual(makedirs.call_count, 0)
        patcher1.stop()
        patcher2.stop()

    def test_FilenameGetterWithDirectoryInitialization_DirectoryDoesNotExist(self):
        settings = self._FilenameGetterSettings(file_name='/foo/kratos.h5')
        patcher1 = patch('os.path.exists', autospec=True)
        patcher2 = patch('os.makedirs', autospec=True)
        pathexists = patcher1.start()
        makedirs = patcher2.start()
        pathexists.return_value = False
        obj = file_io._FilenameGetterWithDirectoryInitialization(settings)
        obj.Get(_SurrogateModelPart())
        pathexists.assert_called_once_with('/foo')
        makedirs.assert_called_once_with('/foo')
        patcher1.stop()
        patcher2.stop()

    def test_Create_Settings(self):
        settings = ParametersWrapper()
        file_io.Create(settings)
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
        io = file_io.Create(settings)
        self.assertIsInstance(io, file_io._HDF5SerialFileIO)
        self.assertTrue(hasattr(io, 'filename_getter'))
        self.assertEqual(io.file_access_mode, 'exclusive')
        if os.name == 'posix':
            self.assertEqual(io.file_driver, 'sec2')
        elif os.name == 'nt':
            self.assertEqual(io.file_driver, 'windows')
        self.assertEqual(io.echo_level, 0)


class TestOperations(KratosUnittest.TestCase):

    def test_CreateNonExistingOperation(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'abcdefg'
        with self.assertRaisesRegex(ValueError, r'"operation_type" has invalid value "abcdefg"'):
            operations.Create(settings)

    def test_Prefix_Literal(self):
        model_part = _SurrogateModelPart()
        prefix = operations.model_part.Prefix('/ModelData', model_part)
        self.assertEqual(prefix, '/ModelData')

    def test_Prefix_NonTerminalTime(self):
        model_part = _SurrogateModelPart()
        prefix = operations.model_part.Prefix('/ModelData-<time>', model_part)
        self.assertEqual(prefix, '/ModelData-1.23456789')

    def test_Prefix_FormattedNonTerminalTime(self):
        model_part = _SurrogateModelPart()
        prefix = operations.model_part.Prefix(
            '/ModelData-<time>', model_part, '0.2f')
        self.assertEqual(prefix, '/ModelData-1.23')

    def test_Prefix_NonTerminalIdentifier(self):
        model_part = _SurrogateModelPart()
        prefix = operations.model_part.Prefix(
            '/<model_part_name>-<time>', model_part)
        self.assertEqual(prefix, '/model_part-1.23456789')

    def test_VariableIO_Settings(self):
        settings1 = ParametersWrapper()
        variable_io = operations.VariableIO(settings1)
        settings2 = variable_io.GetSettings(_SurrogateModelPart())
        for settings in [settings1, settings2]:
            self.assertTrue(settings.Has('prefix'))
            self.assertTrue(settings.Has('list_of_variables'))
            self.assertEqual(settings['prefix'], '/ResultsData')
            self.assertTrue(settings['list_of_variables'].IsArray())
            self.assertEqual(len(settings['list_of_variables']), 0)

    def test_VariableIO_GetSettingsWithNonTerminalPrefix(self):
        input_settings = ParametersWrapper('''
            {
                "prefix": "/ModelData/<model_part_name>/<time>",
                "time_format": "0.2f"
            }
            ''')
        variable_io = operations.VariableIO(input_settings)
        settings = variable_io.GetSettings(_SurrogateModelPart())
        self.assertEqual(settings['prefix'], '/ModelData/model_part/1.23')

    def test_ModelPartOutput(self):
        settings = ParametersWrapper()
        model_part_output = operations.Create(settings)
        self.assertTrue(settings.Has('operation_type'))
        self.assertTrue(settings.Has('prefix'))
        self.assertEqual(settings['operation_type'], 'model_part_output')
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ModelPartIO', autospec=True) as p:
            model_part_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            model_part_output(model_part, hdf5_file)
            p.assert_called_once_with(hdf5_file, '/ModelData')
            model_part_io.WriteModelPart.assert_called_once_with(model_part)

    def test_ModelPartOutput_NonTerminalPrefix(self):
        settings = ParametersWrapper('''
            {
                "operation_type": "model_part_output",
                "prefix": "/ModelData/<model_part_name>/<time>",
                "time_format": "0.2f"
            }
            ''')
        model_part_output = operations.Create(settings)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ModelPartIO', autospec=True) as p:
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            model_part_output(model_part, hdf5_file)
            args, _ = p.call_args
            self.assertEqual(args[1], '/ModelData/model_part/1.23')

    def test_NodalSolutionStepDataOutput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'nodal_solution_step_data_output'
        nodal_solution_step_data_output = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(
            nodal_solution_step_data_output, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepDataIO', autospec=True) as p:
            nodal_solution_step_data_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            nodal_solution_step_data_output(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_solution_step_data_io.WriteNodalResults.call_count, 1)

    def test_NodalSolutionStepDataInput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'nodal_solution_step_data_input'
        nodal_solution_step_data_input = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(
            nodal_solution_step_data_input, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepDataIO', autospec=True) as p:
            nodal_solution_step_data_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            nodal_solution_step_data_input(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_solution_step_data_io.ReadNodalResults.call_count, 1)

    def test_NodalDataValueOutput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'nodal_data_value_output'
        nodal_data_value_output = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(nodal_data_value_output, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalDataValueIO', autospec=True) as p:
            nodal_data_value_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            nodal_data_value_output(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_data_value_io.WriteNodalResults.call_count, 1)

    def test_NodalFlagValueOutput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'nodal_flag_value_output'
        nodal_flag_value_output = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(nodal_flag_value_output, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalFlagValueIO', autospec=True) as p:
            nodal_flag_value_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            nodal_flag_value_output(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_flag_value_io.WriteNodalFlags.call_count, 1)


    def test_NodalDataValueInput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'nodal_data_value_input'
        nodal_data_value_input = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(nodal_data_value_input, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalDataValueIO', autospec=True) as p:
            nodal_data_value_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            nodal_data_value_input(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_data_value_io.ReadNodalResults.call_count, 1)

    def test_NodalFlagValueInput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'nodal_flag_value_input'
        nodal_flag_value_input = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(nodal_flag_value_input, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalFlagValueIO', autospec=True) as p:
            nodal_flag_value_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            nodal_flag_value_input(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_flag_value_io.ReadNodalFlags.call_count, 1)

    def test_ElementDataValueOutput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'element_data_value_output'
        element_data_value_output = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(element_data_value_output, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementDataValueIO', autospec=True) as p:
            element_data_value_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            element_data_value_output(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(element_data_value_io.WriteElementResults.call_count, 1)

    def test_ElementFlagValueOutput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'element_flag_value_output'
        element_flag_value_output = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(element_flag_value_output, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementFlagValueIO', autospec=True) as p:
            element_flag_value_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            element_flag_value_output(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(element_flag_value_io.WriteElementFlags.call_count, 1)

    def test_ElementDataValueInput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'element_data_value_input'
        element_data_value_input = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(element_data_value_input, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementDataValueIO', autospec=True) as p:
            element_data_value_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            element_data_value_input(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(element_data_value_io.ReadElementResults.call_count, 1)

    def test_ElementFlagValueInput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'element_flag_value_input'
        element_flag_value_input = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(element_flag_value_input, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementFlagValueIO', autospec=True) as p:
            element_flag_value_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            element_flag_value_input(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(element_flag_value_io.ReadElementFlags.call_count, 1)

    def test_ConditionDataValueOutput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'condition_data_value_output'
        condition_data_value_output = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(condition_data_value_output, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionDataValueIO', autospec=True) as p:
            condition_data_value_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            condition_data_value_output(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(condition_data_value_io.WriteConditionResults.call_count, 1)

    def test_ConditionFlagValueOutput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'condition_flag_value_output'
        condition_flag_value_output = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(condition_flag_value_output, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionFlagValueIO', autospec=True) as p:
            condition_flag_value_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            condition_flag_value_output(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(condition_flag_value_io.WriteConditionFlags.call_count, 1)

    def test_ConditionDataValueInput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'condition_data_value_input'
        condition_data_value_input = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(condition_data_value_input, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionDataValueIO', autospec=True) as p:
            condition_data_value_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            condition_data_value_input(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(condition_data_value_io.ReadConditionResults.call_count, 1)

    def test_ConditionFlagValueInput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'condition_flag_value_input'
        condition_flag_value_input = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(condition_flag_value_input, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionFlagValueIO', autospec=True) as p:
            condition_flag_value_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            condition_flag_value_input(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(condition_flag_value_io.ReadConditionFlags.call_count, 1)

    def test_PrimalBossakOutput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'primal_bossak_output'
        primal_bossak_output = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(primal_bossak_output, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepBossakIO', autospec=True) as p:
            nodal_solution_step_bossak_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            primal_bossak_output(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_solution_step_bossak_io.WriteNodalResults.call_count, 1)

    def test_PrimalBossakInput(self):
        settings = ParametersWrapper()
        settings['operation_type'] = 'primal_bossak_input'
        primal_bossak_input = operations.Create(settings)
        self.assertTrue(settings.Has('prefix'))
        self.assertTrue(settings.Has('list_of_variables'))
        self.assertTrue(settings['list_of_variables'].IsArray())
        self.assertIsInstance(primal_bossak_input, operations.VariableIO)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepBossakIO', autospec=True) as p:
            nodal_solution_step_bossak_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            primal_bossak_input(model_part, hdf5_file)
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_solution_step_bossak_io.ReadNodalResults.call_count, 1)


class TestControllers(KratosUnittest.TestCase):

    def test_CreateNonExistingController(self):
        settings = ParametersWrapper()
        settings['controller_type'] = 'abcdefg'
        with self.assertRaisesRegex(ValueError, r'"controller_type" has invalid value "abcdefg"'):
            controllers.Create(MagicMock(), MagicMock(), settings)

    def test_DefaultController(self):
        model_part = _SurrogateModelPart()
        io_settings = ParametersWrapper()
        controller_settings = ParametersWrapper()
        with patch('KratosMultiphysics.HDF5Application.core.file_io.KratosHDF5.HDF5FileSerial', autospec=True):
            io = file_io.Create(io_settings)
            controller = controllers.Create(
                model_part, io, controller_settings)
            self.assertTrue(controller_settings.Has('controller_type'))
            self.assertEqual(
                controller_settings['controller_type'], 'default_controller')
            self.assertEqual(controller.model_part, model_part)
            self.assertEqual(controller.io, io)
            operation = MagicMock(spec=operations.ModelPartOutput)
            controller()
            controller.Add(operation)
            for i in range(10):
                controller()
            self.assertEqual(operation.call_count, 10)

    def test_TemporalController_CreateWithDefaults(self):
        model_part = _SurrogateModelPart()
        io = file_io.Create(ParametersWrapper())
        controller_settings = ParametersWrapper()
        controller_settings['controller_type'] = 'temporal_controller'
        controller = controllers.Create(model_part, io, controller_settings)
        self.assertEqual(controller.model_part, model_part)
        self.assertEqual(controller.io, io)
        self.assertEqual(controller.time_frequency, 1.0)
        self.assertEqual(controller.step_frequency, 1)
        self.assertEqual(controller.current_time, 0.0)
        self.assertEqual(controller.current_step, 0)

    def test_TemporalController_CreateWithParameters(self):
        model_part = _SurrogateModelPart()
        io = file_io.Create(ParametersWrapper())
        controller_settings = ParametersWrapper()
        controller_settings['controller_type'] = 'temporal_controller'
        controller_settings['time_frequency'] = 2.0
        controller_settings['step_frequency'] = 3
        controller = controllers.Create(model_part, io, controller_settings)
        self.assertEqual(controller.time_frequency, 2.0)
        self.assertEqual(controller.step_frequency, 3)

    def test_TemporalController_StepFrequency(self):
        model_part = _SurrogateModelPart()
        controller_settings = ParametersWrapper()
        controller_settings['step_frequency'] = 2
        controller_settings['controller_type'] = 'temporal_controller'
        with patch('KratosMultiphysics.HDF5Application.core.file_io._HDF5SerialFileIO', autospec=True):
            io = file_io.Create(ParametersWrapper())
            controller = controllers.Create(
                model_part, io, controller_settings)
            for i in range(10):
                controller()
            io.Get.assert_called_with(model_part)
            self.assertEqual(io.Get.call_count, 5)

    def test_TemporalController_TimeFrequency(self):
        model_part = _SurrogateModelPart()
        controller_settings = ParametersWrapper()
        controller_settings['step_frequency'] = 100
        controller_settings['time_frequency'] = 0.5
        controller_settings['controller_type'] = 'temporal_controller'
        with patch('KratosMultiphysics.HDF5Application.core.file_io._HDF5SerialFileIO', autospec=True):
            io = file_io.Create(ParametersWrapper())
            controller = controllers.Create(
                model_part, io, controller_settings)
            for i in range(10):
                controller()
            io.Get.assert_called_with(model_part)
            self.assertEqual(io.Get.call_count, 2)

    def test_TemporalController_NearlyTheSameTimeFrequency(self):
        model_part = _SurrogateModelPart()
        controller_settings = ParametersWrapper()
        controller_settings['step_frequency'] = 100
        controller_settings['time_frequency'] = 0.2000001
        controller_settings['controller_type'] = 'temporal_controller'
        with patch('KratosMultiphysics.HDF5Application.core.file_io._HDF5SerialFileIO', autospec=True):
            io = file_io.Create(ParametersWrapper())
            controller = controllers.Create(
                model_part, io, controller_settings)
            for _ in range(10):
                controller()
            io.Get.assert_called_with(model_part)
            self.assertEqual(io.Get.call_count, 5)

    def test_TemporalController_OperationCall(self):
        model_part = _SurrogateModelPart()
        controller_settings = ParametersWrapper()
        controller_settings['controller_type'] = 'temporal_controller'
        io = file_io.Create(ParametersWrapper())
        operation = MagicMock(spec=operations.ModelPartOutput)
        controller = controllers.Create(
            model_part, io, controller_settings)
        controller.Add(operation)
        with patch('KratosMultiphysics.HDF5Application.core.file_io.KratosHDF5.HDF5FileSerial', autospec=True):
            for _ in range(10):
                controller()
            self.assertEqual(operation.call_count, 10)


class TestFactory(KratosUnittest.TestCase):

    def test_NonArraySettings(self):
        model = KratosMultiphysics.Model()
        settings = KratosMultiphysics.Parameters()
        with self.assertRaisesRegex(ValueError, r'Expected settings as an array'):
            core.Factory(settings, model)

    def test_EmptyArraySettings(self):
        model = KratosMultiphysics.Model()
        settings = KratosMultiphysics.Parameters('''
            {
                "list_of_controllers" : []
            }
            ''')
        settings = ParametersWrapper(settings)
        with self.assertRaisesRegex(RuntimeError, '"PLEASE_SPECIFY_MODEL_PART_NAME" was not found'):
            core.Factory(settings['list_of_controllers'], model)

    def test_DefaultSettings(self):
        model = KratosMultiphysics.Model()
        model.CreateModelPart('test')
        parent_settings = KratosMultiphysics.Parameters('''
            {
                "list_of_controllers" : [
                    {
                        "model_part_name" : "test"
                    }
                ]
            }
            ''')
        parent_settings = ParametersWrapper(parent_settings)
        core.Factory(parent_settings['list_of_controllers'], model)
        settings = parent_settings['list_of_controllers'][0]
        self.assertTrue(settings.Has('model_part_name'))
        self.assertTrue(settings.Has('process_step'))
        self.assertTrue(settings.Has('controller_settings'))
        self.assertTrue(settings.Has('io_settings'))
        self.assertTrue(settings.Has('list_of_operations'))
        self.assertTrue(settings['list_of_operations'].IsArray())
        self.assertEqual(settings['list_of_operations'].size(), 1)
        self.assertTrue(settings['list_of_operations']
                        [0].Has('operation_type'))
        self.assertTrue(settings['list_of_operations'][0].Has('prefix'))

    def test_DefaultProcess(self):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart('test')
        parent_settings = KratosMultiphysics.Parameters('''
            {
                "list_of_controllers" : [
                    {
                        "model_part_name" : "test"
                    }
                ]
            }
            ''')
        parent_settings = ParametersWrapper(parent_settings)
        process = core.Factory(
            parent_settings['list_of_controllers'], model)
        patcher1 = patch(
            'KratosMultiphysics.HDF5Application.core.file_io.KratosHDF5.HDF5FileSerial', autospec=True)
        patcher2 = patch(
            'KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ModelPartIO', autospec=True)
        patcher1.start()
        MockHDF5ModelPartIO = patcher2.start()
        process.ExecuteInitialize()
        model_part_io = MockHDF5ModelPartIO.return_value
        model_part_io.WriteModelPart.assert_called_once_with(model_part)
        patcher1.stop()
        patcher2.stop()


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
