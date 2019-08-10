import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application import core
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
        settings = KratosMultiphysics.Parameters()
        settings.AddEmptyValue('file_name')
        if 'file_name' in kwargs:
            settings['file_name'].SetString(kwargs['file_name'])
        else:
            settings['file_name'].SetString('kratos.h5')
        if 'time_format' in kwargs:
            settings.AddEmptyValue('time_format').SetString(
                kwargs['time_format'])
        return settings

    def test_FileIO_FileSettings(self):
        io = file_io._FileIO()
        self._BuildTestFileIOObject(io)
        settings = io._FileSettings('kratos.h5')
        self.assertIsInstance(settings, KratosMultiphysics.Parameters)
        self.assertEqual(settings['file_name'].GetString(), 'kratos.h5')
        self.assertEqual(settings['file_access_mode'].GetString(), 'exclusive')
        self.assertEqual(settings['file_driver'].GetString(), 'core')
        self.assertEqual(settings['echo_level'].GetInt(), 0)

    def test_HDF5SerialFileIO_Creation(self):
        io = file_io._HDF5SerialFileIO()
        self._BuildTestFileIOObject(io)
        obj = io.Get('kratos.h5')
        self.assertIsInstance(obj, KratosHDF5.HDF5FileSerial)

    def test_SetDefaults(self):
        settings = KratosMultiphysics.Parameters()
        file_io._SetDefaults(settings)
        self.assertTrue(settings.Has('io_type'))
        self.assertTrue(settings.Has('file_name'))
        self.assertFalse(settings.Has('time_format'))
        self.assertTrue(settings.Has('file_access_mode'))
        self.assertTrue(settings.Has('file_driver'))
        self.assertTrue(settings.Has('echo_level'))
        self.assertEqual(
            settings['io_type'].GetString(), 'serial_hdf5_file_io')
        self.assertEqual(settings['file_name'].GetString(), 'kratos')
        self.assertEqual(settings['file_access_mode'].GetString(), 'exclusive')
        if os.name == 'posix':
            self.assertEqual(settings['file_driver'].GetString(), 'sec2')
        elif os.name == 'nt':
            self.assertEqual(settings['file_driver'].GetString(), 'windows')
        self.assertEqual(settings['echo_level'].GetInt(), 0)

    def test_SetDefaults_NonTerminalTime(self):
        settings = KratosMultiphysics.Parameters()
        settings.AddEmptyValue('file_name').SetString('kratos-<time>.h5')
        file_io._SetDefaults(settings)
        self.assertTrue(settings.Has('time_format'))
        self.assertEqual(settings['time_format'].GetString(), '0.4f')

    def test_SetDefaults_ParallelIO(self):
        settings = KratosMultiphysics.Parameters()
        settings.AddEmptyValue('io_type').SetString('parallel_hdf5_file_io')
        file_io._SetDefaults(settings)
        self.assertEqual(settings['file_driver'].GetString(), 'mpio')

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
        settings = KratosMultiphysics.Parameters()
        file_io.Create(settings)
        self.assertTrue(settings.Has('io_type'))
        self.assertTrue(settings.Has('file_name'))
        self.assertTrue(settings.Has('file_access_mode'))
        self.assertTrue(settings.Has('file_driver'))
        self.assertTrue(settings.Has('echo_level'))
        self.assertEqual(
            settings['io_type'].GetString(), 'serial_hdf5_file_io')
        self.assertEqual(settings['file_access_mode'].GetString(), 'exclusive')
        if os.name == 'posix':
            self.assertEqual(settings['file_driver'].GetString(), 'sec2')
        elif os.name == 'nt':
            self.assertEqual(settings['file_driver'].GetString(), 'windows')
        self.assertEqual(settings['echo_level'].GetInt(), 0)

    def test_Create_Attributes(self):
        settings = KratosMultiphysics.Parameters()
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
        settings = KratosMultiphysics.Parameters()
        settings.AddEmptyValue('operation_type').SetString('abcdefg')
        with self.assertRaisesRegex(ValueError, r'"operation_type" has invalid value "abcdefg"'):
            operations.Create(settings)

    def test_Prefix_Literal(self):
        model_part = _SurrogateModelPart()
        prefix = operations.model_part._Prefix('/ModelData', model_part)
        self.assertEqual(prefix, '/ModelData')

    def test_Prefix_NonTerminalTime(self):
        model_part = _SurrogateModelPart()
        prefix = operations.model_part._Prefix('/ModelData-<time>', model_part)
        self.assertEqual(prefix, '/ModelData-1.23456789')

    def test_Prefix_FormattedNonTerminalTime(self):
        model_part = _SurrogateModelPart()
        prefix = operations.model_part._Prefix(
            '/ModelData-<time>', model_part, '0.2f')
        self.assertEqual(prefix, '/ModelData-1.23')

    def test_Prefix_NonTerminalIdentifier(self):
        model_part = _SurrogateModelPart()
        prefix = operations.model_part._Prefix(
            '/<model_part_name>-<time>', model_part)
        self.assertEqual(prefix, '/model_part-1.23456789')

    def test_VariableIO_Settings(self):
        settings1 = KratosMultiphysics.Parameters()
        variable_io = operations.VariableIO(settings1)
        settings2 = variable_io.GetSettings(_SurrogateModelPart())
        for settings in [settings1, settings2]:
            self.assertTrue(settings.Has('prefix'))
            self.assertTrue(settings.Has('list_of_variables'))
            self.assertEqual(settings['prefix'].GetString(), '/ResultsData')
            self.assertTrue(settings['list_of_variables'].IsArray())
            self.assertEqual(settings['list_of_variables'].size(), 0)

    def test_VariableIO_GetSettingsWithNonTerminalPrefix(self):
        input_settings = KratosMultiphysics.Parameters('''
            {
                "prefix": "/ModelData/<model_part_name>/<time>",
                "time_format": "0.2f"
            }
            ''')
        variable_io = operations.VariableIO(input_settings)
        settings = variable_io.GetSettings(_SurrogateModelPart())
        self.assertEqual(settings['prefix'].GetString(),
                         '/ModelData/model_part/1.23')

    def test_ModelPartOutput(self):
        settings = KratosMultiphysics.Parameters()
        model_part_output = operations.Create(settings)
        self.assertTrue(settings.Has('operation_type'))
        self.assertTrue(settings.Has('prefix'))
        self.assertEqual(
            settings['operation_type'].GetString(), 'model_part_output')
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ModelPartIO', autospec=True) as p:
            model_part_io = p.return_value
            model_part = _SurrogateModelPart()
            hdf5_file = MagicMock(spec=KratosHDF5.HDF5FileSerial)
            model_part_output(model_part, hdf5_file)
            p.assert_called_once_with(hdf5_file, '/ModelData')
            model_part_io.WriteModelPart.assert_called_once_with(model_part)

    def test_ModelPartOutput_NonTerminalPrefix(self):
        settings = KratosMultiphysics.Parameters('''
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
        settings = KratosMultiphysics.Parameters()
        settings.AddEmptyValue('operation_type').SetString(
            'nodal_solution_step_data_output')
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
        settings = KratosMultiphysics.Parameters()
        settings.AddEmptyValue('operation_type').SetString(
            'nodal_solution_step_data_input')
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
        settings = KratosMultiphysics.Parameters()
        settings.AddEmptyValue('operation_type').SetString(
            'nodal_data_value_output')
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

    def test_NodalDataValueInput(self):
        settings = KratosMultiphysics.Parameters()
        settings.AddEmptyValue('operation_type').SetString(
            'nodal_data_value_input')
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

    def test_ElementDataValueOutput(self):
        settings = KratosMultiphysics.Parameters()
        settings.AddEmptyValue('operation_type').SetString(
            'element_data_value_output')
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

    def test_ElementDataValueInput(self):
        settings = KratosMultiphysics.Parameters()
        settings.AddEmptyValue('operation_type').SetString(
            'element_data_value_input')
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

    def test_PrimalBossakOutput(self):
        settings = KratosMultiphysics.Parameters()
        settings.AddEmptyValue('operation_type').SetString(
            'primal_bossak_output')
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
        settings = KratosMultiphysics.Parameters()
        settings.AddEmptyValue('operation_type').SetString(
            'primal_bossak_input')
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
        settings = KratosMultiphysics.Parameters()
        settings.AddEmptyValue('controller_type').SetString('abcdefg')
        with self.assertRaisesRegex(ValueError, r'"controller_type" has invalid value "abcdefg"'):
            controllers.Create(MagicMock(), MagicMock(), settings)

    def test_DefaultController(self):
        model_part = _SurrogateModelPart()
        io_settings = KratosMultiphysics.Parameters()
        controller_settings = KratosMultiphysics.Parameters()
        with patch('KratosMultiphysics.HDF5Application.core.file_io.KratosHDF5.HDF5FileSerial', autospec=True):
            io = file_io.Create(io_settings)
            controller = controllers.Create(
                model_part, io, controller_settings)
            self.assertTrue(controller_settings.Has('controller_type'))
            self.assertEqual(
                controller_settings['controller_type'].GetString(), 'default_controller')
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
        io = file_io.Create(KratosMultiphysics.Parameters())
        controller_settings = KratosMultiphysics.Parameters()
        controller_settings.AddEmptyValue(
            'controller_type').SetString('temporal_controller')
        controller = controllers.Create(model_part, io, controller_settings)
        self.assertEqual(controller.model_part, model_part)
        self.assertEqual(controller.io, io)
        self.assertEqual(controller.time_frequency, 1.0)
        self.assertEqual(controller.step_frequency, 1)
        self.assertEqual(controller.current_time, 0.0)
        self.assertEqual(controller.current_step, 0)

    def test_TemporalController_CreateWithParameters(self):
        model_part = _SurrogateModelPart()
        io = file_io.Create(KratosMultiphysics.Parameters())
        controller_settings = KratosMultiphysics.Parameters()
        controller_settings.AddEmptyValue(
            'controller_type').SetString('temporal_controller')
        controller_settings.AddEmptyValue('time_frequency').SetDouble(2.0)
        controller_settings.AddEmptyValue('step_frequency').SetInt(3)
        controller = controllers.Create(model_part, io, controller_settings)
        self.assertEqual(controller.time_frequency, 2.0)
        self.assertEqual(controller.step_frequency, 3)

    def test_TemporalController_StepFrequency(self):
        model_part = _SurrogateModelPart()
        controller_settings = KratosMultiphysics.Parameters()
        controller_settings.AddEmptyValue('step_frequency').SetInt(2)
        controller_settings.AddEmptyValue(
            'controller_type').SetString('temporal_controller')
        with patch('KratosMultiphysics.HDF5Application.core.file_io._HDF5SerialFileIO', autospec=True):
            io = file_io.Create(KratosMultiphysics.Parameters())
            controller = controllers.Create(
                model_part, io, controller_settings)
            for i in range(10):
                controller()
            io.Get.assert_called_with(model_part)
            self.assertEqual(io.Get.call_count, 5)

    def test_TemporalController_TimeFrequency(self):
        model_part = _SurrogateModelPart()
        controller_settings = KratosMultiphysics.Parameters()
        controller_settings.AddEmptyValue('step_frequency').SetInt(100)
        controller_settings.AddEmptyValue('time_frequency').SetDouble(0.5)
        controller_settings.AddEmptyValue(
            'controller_type').SetString('temporal_controller')
        with patch('KratosMultiphysics.HDF5Application.core.file_io._HDF5SerialFileIO', autospec=True):
            io = file_io.Create(KratosMultiphysics.Parameters())
            controller = controllers.Create(
                model_part, io, controller_settings)
            for i in range(10):
                controller()
            io.Get.assert_called_with(model_part)
            self.assertEqual(io.Get.call_count, 2)

    def test_TemporalController_NearlyTheSameTimeFrequency(self):
        model_part = _SurrogateModelPart()
        controller_settings = KratosMultiphysics.Parameters()
        controller_settings.AddEmptyValue('step_frequency').SetInt(100)
        controller_settings.AddEmptyValue(
            'time_frequency').SetDouble(0.2000001)
        controller_settings.AddEmptyValue(
            'controller_type').SetString('temporal_controller')
        with patch('KratosMultiphysics.HDF5Application.core.file_io._HDF5SerialFileIO', autospec=True):
            io = file_io.Create(KratosMultiphysics.Parameters())
            controller = controllers.Create(
                model_part, io, controller_settings)
            for _ in range(10):
                controller()
            io.Get.assert_called_with(model_part)
            self.assertEqual(io.Get.call_count, 5)

    def test_TemporalController_OperationCall(self):
        model_part = _SurrogateModelPart()
        controller_settings = KratosMultiphysics.Parameters()
        controller_settings.AddEmptyValue(
            'controller_type').SetString('temporal_controller')
        io = file_io.Create(KratosMultiphysics.Parameters())
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
