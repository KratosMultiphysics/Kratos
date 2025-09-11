# --- Core Imports ---
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting

# --- HDF5 Imports ---
from KratosMultiphysics.HDF5Application.core.pattern import EvaluatePattern
from KratosMultiphysics.HDF5Application.core import controllers
from KratosMultiphysics.HDF5Application.core import operations
from KratosMultiphysics.HDF5Application.core import file_io

# --- STD Imports ---
import os
import pathlib
import typing
import functools
from unittest.mock import call, patch, MagicMock


test_file_path = pathlib.Path(".h5")


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
        settings = KratosMultiphysics.Parameters("""{}""")
        if 'file_name' in kwargs:
            settings.AddString("file_name", kwargs['file_name'])
        else:
            settings.AddString("file_name", "kratos.h5")
        if 'time_format' in kwargs:
            settings.AddString("time_format", kwargs['time_format'])
        return settings

    def test_GetFileName(self):
        self.assertEqual(file_io.GetFileName(_SurrogateModelPart()[1], "<model_part_name>-<time>", "0.4f"), "model_part-1.2346.h5")
        self.assertEqual(file_io.GetFileName(_SurrogateModelPart()[1], "<model_part_name>-<time>.h5", "0.4f"), "model_part-1.2346.h5")

    def test_KeepMostRecentFilesWithDelete(self):
        with patch('KratosMultiphysics.HDF5Application.core.file_io.kratos_utils.DeleteFileIfExisting', autospec=True) as p1:
            with patch("KratosMultiphysics.HDF5Application.core.file_io.GetMatchingEntities", autospec=True) as p2:
                p2.return_value = ["file1", "file2", "file3", "file4", "file5", "file6", "file7", "file8"]
                file_io.KeepMostRecentFiles("file<step>", 4, 2)
                self.assertEqual(p1.call_args_list, [call('file3'), call('file4'), call('file5'), call('file6')])

    def test_KeepMostRecentFilesWithOutDelete(self):
        with patch('KratosMultiphysics.HDF5Application.core.file_io.kratos_utils.DeleteFileIfExisting', autospec=True) as p1:
            with patch("KratosMultiphysics.HDF5Application.core.file_io.GetMatchingEntities", autospec=True) as p2:
                p2.return_value = ["file1", "file2", "file3", "file4", "file5", "file6", "file7", "file8"]
                file_io.KeepMostRecentFiles("file<step>", 8, 2)
                self.assertEqual(p1.call_args_list, [])

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
        with self.assertRaises(RuntimeError):
            operations.CreateAggregatedOperation(model, settings)

    @ensure_no_file(test_file_path)
    def test_Prefix_Literal(self):
        _, model_part = _SurrogateModelPart()
        prefix = EvaluatePattern('/ModelData', model_part)
        self.assertEqual(prefix, '/ModelData')

    @ensure_no_file(test_file_path)
    def test_Prefix_NonTerminalTime(self):
        _, model_part = _SurrogateModelPart()
        prefix = EvaluatePattern('/ModelData-<time>', model_part)
        self.assertEqual(prefix, '/ModelData-1.23456789')

    @ensure_no_file(test_file_path)
    def test_Prefix_FormattedNonTerminalTime(self):
        _, model_part = _SurrogateModelPart()
        prefix = EvaluatePattern(
            '/ModelData-<time>', model_part, '0.2f')
        self.assertEqual(prefix, '/ModelData-1.23')

    @ensure_no_file(test_file_path)
    def test_Prefix_NonTerminalIdentifier(self):
        _, model_part = _SurrogateModelPart()
        prefix = EvaluatePattern(
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
        model_part_output = operations.CreateAggregatedOperation(model, settings)
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
        model_part_output = operations.CreateAggregatedOperation(model, settings)
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ModelPartIO', autospec=True) as p:
            model_part_output.Execute()
            args, _ = p.call_args
            current_settings = KratosMultiphysics.Parameters("""{
                "prefix": "/ModelData/model_part/1.23",
                "time_format": "0.2f",
                "custom_attributes":{}
            }""")
            self.assertTrue(current_settings.IsEquivalentTo(args[0]))
            self.assertEqual(args[1].GetFileName(), ".h5")

    @ensure_no_file(test_file_path)
    def test_NodalSolutionStepDataOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "nodal_solution_step_data_output"
            }]
        }""")
        nodal_solution_step_data_output = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepDataIO', autospec=True) as p:
            nodal_solution_step_data_io = p.return_value
            nodal_solution_step_data_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_solution_step_data_io.Write.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_NodalSolutionStepDataInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "nodal_solution_step_data_input"
            }]
        }""")
        nodal_solution_step_data_input = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepDataIO', autospec=True) as p:
            nodal_solution_step_data_io = p.return_value
            nodal_solution_step_data_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_solution_step_data_io.Read.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_NodalDataValueOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "nodal_data_value_output"
            }]
        }""")
        nodal_data_value_output = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalDataValueIO', autospec=True) as p:
            nodal_data_value_io = p.return_value
            nodal_data_value_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_data_value_io.Write.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_NodalFlagValueOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "nodal_flag_value_output"
            }]
        }""")
        nodal_flag_value_output = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalFlagValueIO', autospec=True) as p:
            nodal_flag_value_io = p.return_value
            nodal_flag_value_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_flag_value_io.Write.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_NodalDataValueInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "nodal_data_value_input"
            }]
        }""")
        nodal_data_value_input = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalDataValueIO', autospec=True) as p:
            nodal_data_value_io = p.return_value
            nodal_data_value_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_data_value_io.Read.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_NodalFlagValueInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "nodal_flag_value_input"
            }]
        }""")
        nodal_flag_value_input = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalFlagValueIO', autospec=True) as p:
            nodal_flag_value_io = p.return_value
            nodal_flag_value_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_flag_value_io.Read.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ElementDataValueOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "element_data_value_output"
            }]
        }""")
        element_data_value_output = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementDataValueIO', autospec=True) as p:
            element_data_value_io = p.return_value
            element_data_value_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(element_data_value_io.Write.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ElementFlagValueOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "element_flag_value_output"
            }]
        }""")
        element_flag_value_output = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementFlagValueIO', autospec=True) as p:
            element_flag_value_io = p.return_value
            element_flag_value_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(element_flag_value_io.Write.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ElementDataValueInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "element_data_value_input"
            }]
        }""")
        element_data_value_input = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementDataValueIO', autospec=True) as p:
            element_data_value_io = p.return_value
            element_data_value_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(element_data_value_io.Read.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ElementFlagValueInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "element_flag_value_input"
            }]
        }""")
        element_flag_value_input = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ElementFlagValueIO', autospec=True) as p:
            element_flag_value_io = p.return_value
            element_flag_value_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(element_flag_value_io.Read.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ConditionDataValueOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "condition_data_value_output"
            }]
        }""")
        condition_data_value_output = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionDataValueIO', autospec=True) as p:
            condition_data_value_io = p.return_value
            condition_data_value_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(condition_data_value_io.Write.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ConditionFlagValueOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "condition_flag_value_output"
            }]
        }""")
        condition_flag_value_output = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionFlagValueIO', autospec=True) as p:
            condition_flag_value_io = p.return_value
            condition_flag_value_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(condition_flag_value_io.Write.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ConditionDataValueInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "condition_data_value_input"
            }]
        }""")
        condition_data_value_input = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionDataValueIO', autospec=True) as p:
            condition_data_value_io = p.return_value
            condition_data_value_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(condition_data_value_io.Read.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_ConditionFlagValueInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "condition_flag_value_input"
            }]
        }""")
        condition_flag_value_input = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5ConditionFlagValueIO', autospec=True) as p:
            condition_flag_value_io = p.return_value
            condition_flag_value_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(condition_flag_value_io.Read.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_PrimalBossakOutput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "primal_bossak_output"
            }]
        }""")
        primal_bossak_output = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepBossakIO', autospec=True) as p:
            nodal_solution_step_bossak_io = p.return_value
            primal_bossak_output.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_solution_step_bossak_io.Write.call_count, 1)

    @ensure_no_file(test_file_path)
    def test_PrimalBossakInput(self):
        model, _ = _SurrogateModelPart()
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "model_part",
            "list_of_operations" : [{
                "operation_type" : "primal_bossak_input"
            }]
        }""")
        primal_bossak_input = operations.CreateAggregatedOperation(model, settings)
        operation_settings = settings["list_of_operations"][0]
        self.assertTrue(operation_settings.Has('prefix'))
        self.assertTrue(operation_settings.Has('list_of_variables'))
        self.assertTrue(operation_settings['list_of_variables'].IsArray())
        with patch('KratosMultiphysics.HDF5Application.core.operations.KratosHDF5.HDF5NodalSolutionStepBossakIO', autospec=True) as p:
            nodal_solution_step_bossak_io = p.return_value
            primal_bossak_input.Execute()
            self.assertEqual(p.call_count, 1)
            self.assertEqual(nodal_solution_step_bossak_io.Read.call_count, 1)


class TestControllers(KratosUnittest.TestCase):

    def test_CreateNonExistingController(self):
        settings = KratosMultiphysics.Parameters("""{"controller_type": "abcdefg"}""")
        with self.assertRaises(ValueError):
            controllers.Factory(MagicMock(), settings)

    def test_DefaultController(self):
        controller_settings = KratosMultiphysics.Parameters("""{"controller_type": "default_controller"}""")
        controller = controllers.Factory(KratosMultiphysics.Model(), controller_settings)
        for _ in range(100):
            self.assertTrue(controller.Evaluate())
            controller.Update()

    def test_SingleTimeControllerStep(self):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("test")
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 0

        controller_settings = KratosMultiphysics.Parameters("""
            {
                "controller_type"    : "single_time_controller",
                "model_part_name"    : "test",
                "output_control_type": "step",
                "output_interval"    : 5
            }""")

        # test without restarting
        controller = controllers.Factory(model, controller_settings.Clone())
        controller.Check()
        for i in range(100):
            model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
            if i == 4:
                self.assertTrue(controller.Evaluate())
            else:
                self.assertFalse(controller.Evaluate())
            controller.Update()


        # test with restarting
        controller = controllers.Factory(model, controller_settings.Clone())
        controller.Check()
        for i in range(100):
            model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
            self.assertFalse(controller.Evaluate())
            controller.Update()

    def test_SingleTimeControllerTime(self):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("test")
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0

        controller_settings = KratosMultiphysics.Parameters("""
            {
                "controller_type"    : "single_time_controller",
                "model_part_name"    : "test",
                "output_control_type": "time",
                "output_interval"    : 5
            }""")

        # test without restarting
        controller = controllers.Factory(model, controller_settings.Clone())
        controller.Check()
        for i in range(100):
            model_part.ProcessInfo[KratosMultiphysics.TIME] += 1
            if i == 4:
                self.assertTrue(controller.Evaluate())
            else:
                self.assertFalse(controller.Evaluate())
            controller.Update()


        # test with restarting
        controller = controllers.Factory(model, controller_settings.Clone())
        controller.Check()
        for i in range(100):
            model_part.ProcessInfo[KratosMultiphysics.TIME] += 1
            self.assertFalse(controller.Evaluate())
            controller.Update()

if __name__ == "__main__":
    KratosUnittest.main()
