# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class IgaTestFactory(KratosUnittest.TestCase):
    def setUp(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):

            # Reading the ProjectParameters
            with open(self.file_name + "_project_parameters.json",'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            # To avoid many prints
            if ProjectParameters["problem_data"]["echo_level"].GetInt() == 0:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
            else:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.INFO)

            # Creating the test
            model = KratosMultiphysics.Model()
            self.test = StructuralMechanicsAnalysis(model, ProjectParameters)
            self.test.Initialize()

    def test_execution(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.test.RunSolutionLoop()

    def tearDown(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.test.Finalize()

class SinglePatchTest(IgaTestFactory):
    file_name = "single_patch_test/single_patch"

# Kirchoff-Love Shell
class Shell3pLinearBeamThick(IgaTestFactory):
    file_name = "linear_beam_thick_p4_nCP5/shell_3p"
class Shell3pNonLinearBeamThick(IgaTestFactory):
    file_name = "nonlinear_beam_thick_p2_nCP22/shell_3p"
class Shell3pNonLinearBeamThickSD(IgaTestFactory):
    file_name = "nonlinear_beam_thick_sd_p3_nCP83/shell_3p"
class Shell3pLinearScordelis(IgaTestFactory):
    file_name = "linear_scordelis_p4_nCP10/shell_3p"

if __name__ == '__main__':
    KratosUnittest.main()
