import KratosMultiphysics
try:
  from KratosMultiphysics.StructuralMechanicsApplication.adaptative_remeshing_structural_mechanics_analysis import AdaptativeRemeshingStructuralMechanicsAnalysis
except ImportError as e:
    pass
import os

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import KratosUtilities
from KratosMultiphysics.kratos_utilities import DeleteFilesEndingWith

# TODO: Add fluid counter part
class StructuralMechanicsRemeshingTest(KratosUnittest.TestCase):

    def setUp(self):
        self.work_folder = "mmg_lagrangian_test"

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testTwoDDynamicBeamTest(self):
        self.file_name = "beam2D_test"
        self.__test_execution()

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testThreeDDynamicBeamTest(self):
        self.file_name = "beam3D_test"
        self.__test_execution()

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testTwoDDynamicBeamLineLoadTest(self):
        self.file_name = "beam2D_line_load_test"
        self.__test_execution()

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testThreeDShellTest(self):
        self.file_name = "test_remesh_shell"
        self.__test_execution()

    def tearDown(self):
        DeleteFilesEndingWith(self.work_folder, ".mesh")
        DeleteFilesEndingWith(self.work_folder, ".sol")

    def __test_execution(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):

            # Reading the ProjectParameters
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

            # To avoid many prints
            if project_parameters["problem_data"]["echo_level"].GetInt() == 0:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
            else:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.INFO)

            # Creating and running the test
            model = KratosMultiphysics.Model()
            test = AdaptativeRemeshingStructuralMechanicsAnalysis(model, project_parameters)
            test.Run()

if __name__ == '__main__':
    KratosUnittest.main()