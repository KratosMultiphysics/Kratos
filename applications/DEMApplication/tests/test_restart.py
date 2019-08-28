import KratosMultiphysics as Kratos
from KratosMultiphysics.DEMApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.DEMApplication import DEM_analysis_stage
import os
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, the_type, value, traceback):
        os.chdir(self.currentPath)

class DEMRestartTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        with open("restart_files/ProjectParametersDEM.json", 'r') as parameter_file:
            self.project_parameters_save = Kratos.Parameters(parameter_file.read())

        # To avoid many prints
        # if (self.project_parameters_save["problem_data"]["echo_level"].GetInt() == 0):
        #     Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)

        # # Set common settings
        # self.time_step = 1.5
        # self.start_time = 0.0
        # self.end_time = 5.9

        # self.project_parameters_save["solver_settings"]["time_stepping"]["time_step"].SetDouble(self.time_step)
        # self.project_parameters_save["problem_data"]["start_time"].SetDouble(self.start_time)
        # self.project_parameters_save["problem_data"]["end_time"].SetDouble(self.end_time)

        # Now clone the settings after the common settings are set
        self.project_parameters_load = self.project_parameters_save.Clone()

        # Adding the specific settings (minimal to test the default settings)
        # save_restart_parameters = Kratos.Parameters("""{
        #     "restart_processes" : [
        #         {
        #         "python_module"   : "save_restart_process",
        #         "kratos_module"   : "KratosMultiphysics",
        #         "process_name"    : "SaveRestartProcess",
        #         "Parameters"            : {
        #         }
        #     }]
        #     }""")

        # save_restart_parameters["restart_processes"][0]["Parameters"].AddValue(
        #     "model_part_name",
        #     self.project_parameters_save["solver_settings"]["model_part_name"])


        # self.project_parameters_save.AddValue("output_processes", save_restart_parameters)

        # load_mp_import_settings = self.project_parameters_load["solver_settings"]["model_import_settings"]
        # load_mp_import_settings.AddEmptyValue("restart_load_file_label")
        # load_mp_import_settings["input_type"].SetString("rest")
        # load_mp_import_settings["restart_load_file_label"].SetString("3.0")

        # Correct the path
        # restart_file_path = load_mp_import_settings["input_filename"].GetString()
        # restart_file_path = load_mp_import_settings["input_filename"].SetString(os.path.split(restart_file_path)[1])

    def test_execution(self):
        # Within this location context:

        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model_save = Kratos.Model()
            model_load = Kratos.Model()
            save_analysis = DEM_analysis_stage.DEMAnalysisStage(model_save, self.project_parameters_save)
            save_analysis.mdpas_folder_path = save_analysis.main_path + '/restart_files'
            self.project_parameters_load["solver_settings"]["model_import_settings"]["input_type"].SetString("rest")
            load_analysis = DEM_analysis_stage.DEMAnalysisStage(model_load, self.project_parameters_load)

            def NullFunction():
                pass
            load_analysis.CleanUpOperations = NullFunction
            save_analysis.CleanUpOperations = NullFunction
            save_analysis.Run()
            load_analysis.Run()
            self.AssertEquality(model_save, model_load)

    def AssertEquality(self, model_1, model_2):
        mp_names_1 = {str(k) for k in model_1.GetModelPartNames()}
        mp_names_2 = {str(k) for k in model_2.GetModelPartNames()}
        print('N'*100, mp_names_1)
        print('\n', mp_names_2)
        self.assertEqual(list(mp_names_1), list(mp_names_2))

        for mp_name in model_1.GetModelPartNames():
            mp_save = model_1.GetModelPart(mp_name)
            mp_load = model_2.GetModelPart(mp_name)
            for node_save, node_load in zip(mp_save.Nodes, mp_load.Nodes):
                displacement_save = node_save.GetSolutionStepValue(Kratos.DISPLACEMENT)
                displacement_load = node_load.GetSolutionStepValue(Kratos.DISPLACEMENT)
                for d1, d2 in zip(displacement_save, displacement_load):
                    self.assertAlmostEqual(d1, d2)

class TestRestart(DEMRestartTestFactory):
    file_name = "restart_files/two_balls"


if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRestart]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.main()
