import KratosMultiphysics
try:
  from KratosMultiphysics.StructuralMechanicsApplication.adaptative_remeshing_structural_mechanics_analysis import AdaptativeRemeshingStructuralMechanicsAnalysis
except ImportError as e:
    pass
import os

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

# TODO: Add fluid counter part
class MeshingStructuralMechanicsTestFactory(KratosUnittest.TestCase):
    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def setUp(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):

            # Reading the ProjectParameters
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            self.modify_parameters(ProjectParameters)

            # To avoid many prints
            if ProjectParameters["problem_data"]["echo_level"].GetInt() == 0:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
            else:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.INFO)

            # Creating the test
            model = KratosMultiphysics.Model()
            self.test = AdaptativeRemeshingStructuralMechanicsAnalysis(model, ProjectParameters)
            self.test.Initialize()

    def modify_parameters(self, project_parameters):
        """This function can be used in derived classes to modify existing parameters
        before the execution of the test (e.g. switch to MPI)
        """
        pass

    def test_execution(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.test.RunSolutionLoop()

    def tearDown(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.test.Finalize()

class TwoDDynamicBeamTest(MeshingStructuralMechanicsTestFactory):
    file_name = "mmg_lagrangian_test/beam2D_test"

class TwoDDynamicBeamLineLoadTest(MeshingStructuralMechanicsTestFactory):
    file_name = "mmg_lagrangian_test/beam2D_line_load_test"

class ThreeDShellTest(MeshingStructuralMechanicsTestFactory):
    file_name = "mmg_lagrangian_test/test_remesh_shell"

class ThreeDDynamicBeamTest(MeshingStructuralMechanicsTestFactory):
    file_name = "mmg_lagrangian_test/beam3D_test"
