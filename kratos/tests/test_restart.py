from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
import restart_utility

import os
import sys

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def ReadModelPart(file_path):
    model_part_name = "MainRestart"
    model_part = KratosMultiphysics.ModelPart(model_part_name)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
    model_part_io = KratosMultiphysics.ModelPartIO(file_path)
    model_part_io.ReadModelPart(model_part)
    return model_part

class TestRestart(KratosUnittest.TestCase):

    def setUp(self):
        if (sys.version_info < (3, 2)):
            self.assertRaisesRegex = self.assertRaisesRegexp

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting("test_restart_file.rest")
        kratos_utils.DeleteFileIfExisting("test_restart_file_5.3.rest")

    def _check_modelpart(self, model_part):
        self.assertEqual(model_part.NumberOfSubModelParts(), 2)

        self.assertEqual(model_part.NumberOfTables(), 1)
        table = model_part.GetTable(1)
        self.assertAlmostEqual(table.GetValue(250), 2.5e-6, 12)

        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfNodes(), 6)
        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfElements(), 4)
        self.assertEqual(model_part.NumberOfConditions(), 5)

        self.assertTrue(model_part.GetNode(1).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertTrue(model_part.GetNode(2).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertFalse(model_part.GetNode(3).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertFalse(model_part.GetNode(972).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertTrue(model_part.GetNode(973).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertTrue(model_part.GetNode(974).IsFixed(KratosMultiphysics.DISPLACEMENT_X))

        self.assertTrue(model_part.GetNode(1).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertTrue(model_part.GetNode(2).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertFalse(model_part.GetNode(3).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertFalse(model_part.GetNode(972).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertTrue(model_part.GetNode(973).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertTrue(model_part.GetNode(974).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))

        self.assertTrue(model_part.GetNode(1).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertTrue(model_part.GetNode(2).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertFalse(model_part.GetNode(3).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertFalse(model_part.GetNode(972).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertTrue(model_part.GetNode(973).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertTrue(model_part.GetNode(974).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.1)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.2)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.0)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.0)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),0.0)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),0.0)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),0.000973)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),0.000974)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z),0.0)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VISCOSITY),0.01)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VISCOSITY),0.01)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.VISCOSITY),0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.VISCOSITY),0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.VISCOSITY),0.01)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.VISCOSITY),0.01)

        self.assertTrue(model_part.HasSubModelPart("Inlets"))

        inlets_model_part = model_part.GetSubModelPart("Inlets")

        self.assertEqual(inlets_model_part.NumberOfTables(), 1)
        self.assertEqual(inlets_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlets_model_part.NumberOfNodes(), 3)
        self.assertEqual(inlets_model_part.NumberOfElements(), 1)
        self.assertEqual(inlets_model_part.NumberOfConditions(), 3)
        self.assertEqual(inlets_model_part.NumberOfSubModelParts(), 2)
        self.assertTrue(inlets_model_part.HasSubModelPart("Inlet1"))
        self.assertTrue(inlets_model_part.HasSubModelPart("Inlet2"))

        inlet1_model_part = inlets_model_part.GetSubModelPart("Inlet1")

        self.assertEqual(inlet1_model_part.NumberOfTables(), 0)
        self.assertEqual(inlet1_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlet1_model_part.NumberOfNodes(), 2)
        self.assertEqual(inlet1_model_part.NumberOfElements(), 0)
        self.assertEqual(inlet1_model_part.NumberOfConditions(), 2)
        self.assertEqual(inlet1_model_part.NumberOfSubModelParts(), 0)

        inlet2_model_part = inlets_model_part.GetSubModelPart("Inlet2")

        self.assertEqual(inlet2_model_part.NumberOfTables(), 0)
        self.assertEqual(inlet2_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlet2_model_part.NumberOfNodes(), 0)
        self.assertEqual(inlet2_model_part.NumberOfElements(), 0)
        self.assertEqual(inlet2_model_part.NumberOfConditions(), 2)
        self.assertEqual(inlet2_model_part.NumberOfSubModelParts(), 0)

        self.assertTrue(model_part.HasSubModelPart("Outlet"))

        outlet_model_part = model_part.GetSubModelPart("Outlet")

        self.assertEqual(outlet_model_part.NumberOfTables(), 0)
        self.assertEqual(outlet_model_part.NumberOfProperties(), 1)
        self.assertEqual(outlet_model_part.NumberOfNodes(), 0)
        self.assertEqual(outlet_model_part.NumberOfElements(), 0)
        self.assertEqual(outlet_model_part.NumberOfConditions(), 1)
        self.assertEqual(outlet_model_part.NumberOfSubModelParts(), 0)

    def __execute_restart_save(self, file_name, serializer_flag):
        model_part = ReadModelPart(GetFilePath("test_model_part_io_read"))

        serializer_save = KratosMultiphysics.Serializer(file_name, serializer_flag)
        serializer_save.Save(model_part.Name, model_part)

    def __execute_restart_load(self, file_name, serializer_flag):
        model_part_name = "MainRestart"

        loaded_model_part = KratosMultiphysics.ModelPart(model_part_name)

        serializer_load = KratosMultiphysics.Serializer(file_name, serializer_flag)
        serializer_load.Load(loaded_model_part.Name, loaded_model_part)

        return loaded_model_part

    def __execute_restart_test(self, serializer_flag):
        file_name = "test_restart_file"
        self.__execute_restart_save(file_name, serializer_flag)
        model_part = self.__execute_restart_load(file_name, serializer_flag)

        self._check_modelpart(model_part)

    def __execute_restart_utility_save(self, model_part_name, restart_time):
        model_part = ReadModelPart(GetFilePath("test_model_part_io_read"))

        model_part.ProcessInfo[KratosMultiphysics.TIME] = restart_time # saving is only done if time > 0.0

        restart_parameters = KratosMultiphysics.Parameters("""
        {
            "input_filename"                    : "test_restart_file",
            "restart_save_frequency"            : 0.0,
            "save_restart_files_in_folder"      : false
        }
        """)

        rest_utility = restart_utility.RestartUtility(model_part, restart_parameters)

        rest_utility.SaveRestart()

    def __execute_restart_utility_load(self, model_part_name, restart_time):
        loaded_model_part = KratosMultiphysics.ModelPart(model_part_name)

        restart_parameters = KratosMultiphysics.Parameters("""
        {
            "input_filename"                    : "test_restart_file",
            "restart_load_file_label"           : "",
            "load_restart_files_from_folder"    : false
        }
        """)

        restart_parameters["restart_load_file_label"].SetString(str(restart_time))

        rest_utility = restart_utility.RestartUtility(loaded_model_part, restart_parameters)

        rest_utility.LoadRestart()

        return loaded_model_part


    def test_restart_NOTRACE(self):
        self.__execute_restart_test(KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE)

    def test_restart_TRACE_ERROR(self):
        self.__execute_restart_test(KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ERROR)

    def test_restart_TRACE_ALL(self):
        self.__execute_restart_test(KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL)

    def test_restart_utility(self):
        # Here we only test SERIALIZER_NO_TRACE since the others are tested in the simple tests
        model_part_name = "MainRestart"
        restart_time = 5.3
        self.__execute_restart_utility_save(model_part_name, restart_time)
        loaded_model_part = self.__execute_restart_utility_load(model_part_name, restart_time)

        self._check_modelpart(loaded_model_part)


if __name__ == '__main__':
    KratosUnittest.main()
