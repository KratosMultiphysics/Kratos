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

class MembraneSinglePatchFourPointSailLinearStatic(IgaTestFactory):
    file_name = "membrane_test/Membrane_single_patch_four_point_sail_linear_static"

class MembraneSinglePatchFourPointSailNonLinearStatic(IgaTestFactory):
    file_name = "membrane_test/Membrane_single_patch_four_point_sail_non_linear_static"

class MembraneSinglePatchFourPointSailImplicitDynamic(IgaTestFactory):
    file_name = "membrane_test/Membrane_single_patch_four_point_sail_implicit_dynamic"

class MembraneMultiPatchFourPointSailLinearStatic(IgaTestFactory):
    file_name = "membrane_test/Membrane_multi_patch_four_point_sail_linear_static"

class MembraneMultiPatchFourPointSailNonLinearStatic(IgaTestFactory):
    file_name = "membrane_test/Membrane_multi_patch_four_point_sail_non_linear_static"

class MembraneMultiPatchFourPointSailImplicitDynamic(IgaTestFactory):
    file_name = "membrane_test/Membrane_multi_patch_four_point_sail_implicit_dynamic"

class FormfindingMembraneSinglePatchFourPointSail(IgaTestFactory):
    file_name = "formfinding_test/Formfinding_membrane_single_patch_four_point_sail"

class FormfindingMembraneMultiPatchFourPointSail(IgaTestFactory):
    file_name = "formfinding_test/Formfinding_membrane_multi_patch_four_point_sail"
class ScordelisRoofShell3pTest(IgaTestFactory):
    file_name = "scordelis_roof_shell_3p_test/scordelis_roof_shell_3p"

if __name__ == '__main__':
    KratosUnittest.main()
