from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

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

class SinglePatchFourPointSailLinearStatic(IgaTestFactory):
    file_name = "single_patch_four_point_sail_linear_static_test/single_patch_four_point_sail_linear_static"

class SinglePatchFourPointSailNonLinearStatic(IgaTestFactory):
    file_name = "single_patch_four_point_sail_non_linear_static_test/single_patch_four_point_sail_non_linear_static"

class SinglePatchFourPointSailFormFinding(IgaTestFactory):
    file_name = "single_patch_four_point_sail_form_finding_test/single_patch_four_point_sail_form_finding"

if __name__ == '__main__':
    KratosUnittest.main()
