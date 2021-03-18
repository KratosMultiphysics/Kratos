# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication
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

# 3p Kirchhoff-Love Shell
class ScordelisRoofShell3pTest(IgaTestFactory):
    file_name = "scordelis_roof_shell_3p_test/scordelis_roof_shell_3p"

class LinearBeamShell3pTest(IgaTestFactory):
    file_name = "linear_beam_shell_3p_test/linear_beam_shell_3p"

# Hierarchic 5p Shell
class Shell5pHierarchicLinearThickBeamTest(IgaTestFactory):
    file_name = "linear_beam_thick_p4_nCP5/shell_5p"

class Shell5pHierarchicLinearScordelisTest(IgaTestFactory):
    file_name = "linear_scordelis_p4_nCP10/shell_5p"

class Shell5pHierarchicNonLinearThickBeamTest(IgaTestFactory):
    file_name = "nonlinear_beam_thick_p2_nCP22/shell_5p"

if __name__ == '__main__':
    KratosUnittest.main()
