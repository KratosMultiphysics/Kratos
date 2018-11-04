from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
import restart_utility
import save_restart_process as save_rest_proc

import os
import sys

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def ReadModelPart(file_path, current_model):
    model_part_name = "MainRestart"
    model_part = current_model.CreateModelPart(model_part_name)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
    model_part_io = KratosMultiphysics.ModelPartIO(file_path)
    model_part_io.ReadModelPart(model_part)

    # Manually adding a constraint to check the serialization of constraints in the ModelPart
    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VISCOSITY, model_part)
    c1 = KratosMultiphysics.MasterSlaveConstraint(10)
    model_part.AddMasterSlaveConstraint(c1)

    return model_part

def IsRestartFile(file_name):
    return os.path.isfile(file_name) and file_name.endswith('.rest')

class TestRestart(KratosUnittest.TestCase):

    def setUp(self):
        if (sys.version_info < (3, 2)):
            self.assertRaisesRegex = self.assertRaisesRegexp

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting("test_restart_file.rest")
        kratos_utils.DeleteFileIfExisting("test_restart_file_15.0.rest")
        kratos_utils.DeleteDirectoryIfExisting("MainRestart__restart_files")

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

        self.assertTrue( 10 in model_part.MasterSlaveConstraints )

    def __execute_restart_save(self, file_name, serializer_flag):

        temporary_model = KratosMultiphysics.Model() #this lives only until the end of this function

        model_part = ReadModelPart(GetFilePath("test_model_part_io_read"), temporary_model)

        serializer_save = KratosMultiphysics.Serializer(file_name, serializer_flag)
        serializer_save.Save(model_part.Name, model_part)
        

    def __execute_restart_load(self, current_model, file_name, serializer_flag):
        model_part_name = "MainRestart"

        loaded_model_part = current_model.CreateModelPart(model_part_name)

        serializer_load = KratosMultiphysics.Serializer(file_name, serializer_flag)
        serializer_load.Load(loaded_model_part.Name, loaded_model_part)

        return loaded_model_part

    def __execute_restart_test(self, serializer_flag):
        file_name = "test_restart_file"
        self.__execute_restart_save(file_name, serializer_flag)

        current_model = KratosMultiphysics.Model() #here we create the Model which will live to the end of this function
        model_part = self.__execute_restart_load(current_model, file_name, serializer_flag)
        self._check_modelpart(model_part)

    def __execute_restart_utility_save(self, model_part_name, restart_time):
        #creating a Model which will be destroyed at the end of the save function
        temporary_model = KratosMultiphysics.Model()

        model_part = ReadModelPart(GetFilePath("test_model_part_io_read"), temporary_model)

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.0 # saving is only done if time > 0.0

        restart_parameters = KratosMultiphysics.Parameters("""
        {
            "input_filename"               : "test_restart_file",
            "restart_save_frequency"       : 10.0,
            "save_restart_files_in_folder" : false
        }
        """)

        rest_utility = restart_utility.RestartUtility(model_part, restart_parameters)

        self.assertFalse(rest_utility.IsRestartOutputStep())

        model_part.ProcessInfo[KratosMultiphysics.TIME] = restart_time

        self.assertTrue(rest_utility.IsRestartOutputStep())

        if rest_utility.IsRestartOutputStep():
            rest_utility.SaveRestart()

        del temporary_model ##explicitly deleting to be sure memory is freed

    def __execute_restart_utility_load(self, current_model, model_part_name, restart_time):
        loaded_model_part = current_model.CreateModelPart(model_part_name)

        restart_parameters = KratosMultiphysics.Parameters("""
        {
            "input_filename"                 : "test_restart_file",
            "restart_load_file_label"        : "",
            "load_restart_files_from_folder" : false
        }
        """)

        restart_parameters["restart_load_file_label"].SetString(str(restart_time))

        rest_utility = restart_utility.RestartUtility(loaded_model_part, restart_parameters)

        rest_utility.LoadRestart() #TODO: it would be best to return the loaded_modelpart from this... 

        return loaded_model_part #rest_utility.model_part



    def test_restart_NOTRACE(self):
        self.__execute_restart_test(KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE)

    def test_restart_TRACE_ERROR(self):
        self.__execute_restart_test(KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ERROR)

    def test_restart_TRACE_ALL(self):
        self.__execute_restart_test(KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL)

    def test_restart_utility(self):
        

        # Here we only test SERIALIZER_NO_TRACE since the others are tested in the simple tests
        model_part_name = "MainRestart"
        restart_time = 15.0

        self.__execute_restart_utility_save(model_part_name, restart_time)
        #we create here the model to which we will load
        current_model = KratosMultiphysics.Model()
        loaded_model_part = self.__execute_restart_utility_load(current_model, model_part_name, restart_time)

        self._check_modelpart(loaded_model_part)

    def test_save_restart_process(self):
        model = KratosMultiphysics.Model()
        model_part = ReadModelPart(GetFilePath("test_model_part_io_read"), model)

        # Here "step" is used as control type, since "time" (=> default) is covered in the tests above
        save_restart_process_params = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name"        : "MainRestart",
                "restart_save_frequency" : 2,
                "restart_control_type"   : "step"
            }
        }""")

        model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.0
        model_part.ProcessInfo[KratosMultiphysics.STEP] = 0

        delta_time = 0.35
        end_time = 17.1

        save_restart_process = save_rest_proc.Factory(save_restart_process_params, model)
        save_restart_process.ExecuteInitialize()
        save_restart_process.ExecuteBeforeSolutionLoop()
        while model_part.ProcessInfo[KratosMultiphysics.TIME] < end_time:
            model_part.ProcessInfo[KratosMultiphysics.TIME] += delta_time
            model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
            save_restart_process.ExecuteInitializeSolutionStep()
            if save_restart_process.IsOutputStep():
                save_restart_process.PrintOutput()
            save_restart_process.ExecuteFinalizeSolutionStep()
        save_restart_process.ExecuteFinalize()

        # Checking if the files exist
        base_path = "MainRestart__restart_files"
        base_file_name = os.path.join(base_path, "MainRestart_")
        for i in range(2,50,2):
            self.assertTrue(os.path.isfile(base_file_name + str(i) + ".rest"))

        # Check number of restart-files
        expected_num_files = 24
        num_files = len([name for name in os.listdir(base_path) if IsRestartFile(os.path.join(base_path, name))])
        self.assertEqual(expected_num_files, num_files)

        #deleting the model so to be sure nothing remains in the memory
        del(model)


        # Loading one of the files and checking if the loaded model_part is ok
        loaded_model = KratosMultiphysics.Model()
        file_name = base_file_name + "16"
        loaded_model_part = self.__execute_restart_load(loaded_model, file_name, KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE)

        self._check_modelpart(loaded_model_part)


if __name__ == '__main__':
    KratosUnittest.main()
