# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
import run_cpp_unit_tests
import KratosMultiphysics.kratos_utilities as kratos_utilities

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import Iga test factory tests
from iga_test_factory import SinglePatchTest as SinglePatchTest
# Membrane tests
from iga_test_factory import MembraneSinglePatchFourPointSailLinearStatic as MembraneSinglePatchFourPointSailLinearStatic
from iga_test_factory import MembraneSinglePatchFourPointSailNonLinearStatic as MembraneSinglePatchFourPointSailNonLinearStatic
from iga_test_factory import MembraneSinglePatchFourPointSailImplicitDynamic as MembraneSinglePatchFourPointSailImplicitDynamic
# 3p Shell KL
from iga_test_factory import ScordelisRoofShell3pTest as ScordelisRoofShell3pTest
from iga_test_factory import LinearBeamShell3pTest as LinearBeamShell3pTest
# 5p Shell Hierarchic
from iga_test_factory import Shell5pHierarchicLinearThickBeamTest as TShell5pHierarchicLinearThickBeamTest
from iga_test_factory import Shell5pHierarchicLinearScordelisTest as TShell5pHierarchicLinearScordelisTest
from iga_test_factory import Shell5pHierarchicNonLinearThickBeamTest as TShell5pHierarchicNonLinearThickBeamTest
# 5p Shell Hierarchic
from iga_test_factory import ScordelisRoofShell5pTest as ScordelisRoofShell5pTest
# Weak support tests
from iga_test_factory import SinglePatchRefinedSupportPenaltyTest as SinglePatchRefinedSupportPenaltyTest
from iga_test_factory import SinglePatchRefinedSupportLagrangeTest as SinglePatchRefinedSupportLagrangeTest
from iga_test_factory import SinglePatchRefinedSupportNitscheTest as SinglePatchRefinedSupportNitscheTest

# Modelers tests
from test_modelers import TestModelers as TTestModelers
# Nurbs Geometry tests
from test_nurbs_volume_element import TestNurbsVolumeElement as TTestNurbsVolumeElements

has_linear_solvers_application = kratos_utilities.CheckIfApplicationsAvailable("LinearSolversApplication")

def AssembleTestSuites():
    ''' Populates the test suites to run.
    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"
    Return
    ------
    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([
        # Import test
        SinglePatchTest,
        # Membrane tests
        MembraneSinglePatchFourPointSailLinearStatic,
        MembraneSinglePatchFourPointSailNonLinearStatic,
        # 3p Shell KL
        ScordelisRoofShell3pTest,
        LinearBeamShell3pTest,
        # 5p Shell Director
        #ScordelisRoofShell5pTest,
        TTestModelers,
        TTestNurbsVolumeElements,
        # Weak support tests
        SinglePatchRefinedSupportPenaltyTest,
        SinglePatchRefinedSupportLagrangeTest
        ]))

    if has_linear_solvers_application:
        from KratosMultiphysics import LinearSolversApplication
        if LinearSolversApplication.HasFEAST():
            smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([
                # Weak support Nitsche test
                SinglePatchRefinedSupportNitscheTest
                ]))
        else:
            print("FEAST not available in LinearSolversApplication")

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([
        # Membrane tests
        MembraneSinglePatchFourPointSailImplicitDynamic,
        # 5p Shell Hierarchic
        TShell5pHierarchicLinearThickBeamTest,
        TShell5pHierarchicLinearScordelisTest,
        TShell5pHierarchicNonLinearThickBeamTest
        ]))

    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
