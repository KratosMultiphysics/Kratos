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

# Membrane
class MembraneSinglePatchFourPointSailLinearStatic(IgaTestFactory):
    file_name = "membrane_test/Membrane_single_patch_four_point_sail_linear_static"

class MembraneSinglePatchFourPointSailNonLinearStatic(IgaTestFactory):
    file_name = "membrane_test/Membrane_single_patch_four_point_sail_non_linear_static"

class MembraneSinglePatchFourPointSailImplicitDynamic(IgaTestFactory):
    file_name = "membrane_test/Membrane_single_patch_four_point_sail_implicit_dynamic"

# 3p Kirchhoff-Love Shell
@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class ScordelisRoofShell3pTest(IgaTestFactory):
    file_name = "scordelis_roof_test/scordelis_roof_shell_3p"

class LinearBeamShell3pTest(IgaTestFactory):
    file_name = "linear_beam_shell_3p_test/linear_beam_shell_3p"

class LinearBeamShell3pAdditiveSchwarzTest(IgaTestFactory):
    file_name = "linear_beam_shell_3p_test/linear_beam_shell_3p_additive_schwarz"

# Hierarchic 5p Shell
class Shell5pHierarchicLinearThickBeamTest(IgaTestFactory):
    file_name = "linear_beam_thick_p4_nCP5/shell_5p"

class Shell5pHierarchicLinearScordelisTest(IgaTestFactory):
    file_name = "linear_scordelis_p4_nCP10/shell_5p"

class Shell5pHierarchicNonLinearThickBeamTest(IgaTestFactory):
    file_name = "nonlinear_beam_thick_p2_nCP22/shell_5p"

# 5p Shell Director
class ScordelisRoofShell5pTest(IgaTestFactory):
    file_name = "scordelis_roof_test/scordelis_roof_shell_5p"

# Weak support
@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class SinglePatchRefinedSupportPenaltyTest(IgaTestFactory):
    file_name = "weak_support_tests/single_patch_refined_test/single_patch_refined_test_penalty"

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class SinglePatchRefinedSupportLagrangeTest(IgaTestFactory):
    file_name = "weak_support_tests/single_patch_refined_test/single_patch_refined_test_lagrange"

class SinglePatchRefinedSupportNitscheTest(IgaTestFactory):
    file_name = "weak_support_tests/single_patch_refined_test/single_patch_refined_test_nitsche"

# Coupling C_0
@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class TwoPatchCouplingPenaltyShell3pTest(IgaTestFactory):
    file_name = "coupling_condition_tests/two_patch_test/two_patch_test_penalty_shell_3p"

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class TwoPatchCouplingLagrangeShell3pTest(IgaTestFactory):
    file_name = "coupling_condition_tests/two_patch_test/two_patch_test_lagrange_shell_3p"

class TwoPatchCouplingNitscheShell3pTest(IgaTestFactory):
    file_name = "coupling_condition_tests/two_patch_test/two_patch_test_nitsche_shell_3p"

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class TwoPatchRefinedCouplingPenaltyMembraneTest(IgaTestFactory):
    file_name = "coupling_condition_tests/two_patch_refined_test/two_patch_refined_test_penalty_membrane"

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class TwoPatchRefinedCouplingLagrangeMembraneTest(IgaTestFactory):
    file_name = "coupling_condition_tests/two_patch_refined_test/two_patch_refined_test_lagrange_membrane"

class TwoPatchRefinedCouplingNitscheMembraneTest(IgaTestFactory):
    file_name = "coupling_condition_tests/two_patch_refined_test/two_patch_refined_test_nitsche_membrane"

# Rotation/G_1 Coupling
class TwoPatchCantileverCouplingPenaltyTest(IgaTestFactory):
    file_name = "coupling_condition_tests/two_patch_cantilever_test/two_patch_cantilever_test_penalty"

class TwoPatchCantileverRefinedCouplingPenaltyTest(IgaTestFactory):
    file_name = "coupling_condition_tests/two_patch_cantilever_refined_test/two_patch_cantilever_refined_test_penalty"

if __name__ == '__main__':
    KratosUnittest.main()
