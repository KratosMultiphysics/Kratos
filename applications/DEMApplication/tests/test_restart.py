import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.DEMApplication import DEM_analysis_stage

# Defining a generic Test, that is not actually KratosUnittest
class DEMRestartTestFactory():

    def setUp(self, case_name=""):
        with KratosUnittest.WorkFolderScope(".", __file__):
            with open('restart_files/' + case_name + 'ProjectParametersDEM.json', 'r') as parameter_file:
                self.project_parameters_save = Kratos.Parameters(parameter_file.read())

        # Now clone the settings after the common settings are set
        self.project_parameters_load = self.project_parameters_save.Clone()
        self.project_parameters_load["solver_settings"]["model_import_settings"]["input_type"].SetString("rest")

    def test_execution(self):
        # Within this location context:

        with KratosUnittest.WorkFolderScope(".", __file__):
            model_save = Kratos.Model()
            model_load = Kratos.Model()
            save_analysis = DEM_analysis_stage.DEMAnalysisStage(model_save, self.project_parameters_save)
            save_analysis.mdpas_folder_path = save_analysis.main_path + '/restart_files'
            load_analysis = DEM_analysis_stage.DEMAnalysisStage(model_load, self.project_parameters_load)
            self.project_parameters_load["output_processes"].RemoveValue("restart_processes")
            def NullFunction():
                pass
            load_analysis.CleanUpOperations = NullFunction
            save_analysis.CleanUpOperations = NullFunction
            save_analysis.Run()
            load_analysis.Run()
            self.AssertEquality(model_save, model_load)

    def AssertEquality(self, model_1, model_2):
        mp_names_1 = set(model_1.GetModelPartNames())
        mp_names_2 = set(model_2.GetModelPartNames())

        self.assertEqual(mp_names_1, mp_names_2)

        for mp_name in model_1.GetModelPartNames():
            mp_save = model_1.GetModelPart(mp_name)
            mp_load = model_2.GetModelPart(mp_name)
            for node_save, node_load in zip(mp_save.Nodes, mp_load.Nodes):
                displacement_save = node_save.GetSolutionStepValue(Kratos.DISPLACEMENT)
                displacement_load = node_load.GetSolutionStepValue(Kratos.DISPLACEMENT)
                for d1, d2 in zip(displacement_save, displacement_load):
                    self.assertAlmostEqual(d1, d2)

class TestRestartOneBall(DEMRestartTestFactory, KratosUnittest.TestCase):
    case_name = "one_ball"
    def setUp(self):
        super(TestRestartOneBall, self).setUp(TestRestartOneBall.case_name)

class TestRestartTwoBalls(DEMRestartTestFactory, KratosUnittest.TestCase):
    case_name = "two_balls"
    def setUp(self):
        super(TestRestartTwoBalls, self).setUp(TestRestartTwoBalls.case_name)

class TestRestartBallAndWall(DEMRestartTestFactory, KratosUnittest.TestCase):
    case_name = "ball_and_wall"
    def setUp(self):
        super(TestRestartBallAndWall, self).setUp(TestRestartBallAndWall.case_name)

if __name__ == '__main__':
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)

    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    test_list = [
        TestRestartOneBall,
        TestRestartTwoBalls,
        TestRestartBallAndWall
    ]
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(test_list))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.main()
