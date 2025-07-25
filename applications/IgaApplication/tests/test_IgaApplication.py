# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
import KratosMultiphysics.kratos_utilities as kratos_utilities

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import Iga test factory tests
from iga_test_factory import SinglePatchTest as SinglePatchTest
# Truss tests - python based
from truss_element_tests import TrussElementTests as TTrussElementTests
# Structural Elements test - python based
from test_solid_IGA_element import SolidIGAElementTests as TSolidIGAElementTests
# Sbm Structural Elements test
from sbm_solid_tests import SbmSolidTests as TSbmSolidTests

# Membrane tests
from iga_test_factory import MembraneSinglePatchFourPointSailLinearStatic as MembraneSinglePatchFourPointSailLinearStatic
from iga_test_factory import MembraneSinglePatchFourPointSailNonLinearStatic as MembraneSinglePatchFourPointSailNonLinearStatic
from iga_test_factory import MembraneSinglePatchFourPointSailImplicitDynamic as MembraneSinglePatchFourPointSailImplicitDynamic
# 3p Shell KL - python based
from shell_3p_element_tests import Shell3pElementTests as TShell3pElementTests
# 3p Shell KL
from iga_test_factory import ScordelisRoofShell3pTest as ScordelisRoofShell3pTest
from iga_test_factory import LinearBeamShell3pTest as LinearBeamShell3pTest
# 3p Shell Kl, solved with additive schwarz preconditioner
from iga_test_factory import LinearBeamShell3pAdditiveSchwarzTest as LinearBeamShell3pAdditiveSchwarzTest
# 5p Shell Hierarchic
from iga_test_factory import Shell5pHierarchicLinearThickBeamTest as TShell5pHierarchicLinearThickBeamTest
from iga_test_factory import Shell5pHierarchicLinearScordelisTest as TShell5pHierarchicLinearScordelisTest
from iga_test_factory import Shell5pHierarchicNonLinearThickBeamTest as TShell5pHierarchicNonLinearThickBeamTest
# 5p Shell
from iga_test_factory import ScordelisRoofShell5pTest as ScordelisRoofShell5pTest
# Weak support tests
from iga_test_factory import SinglePatchRefinedSupportPenaltyTest as SinglePatchRefinedSupportPenaltyTest
from iga_test_factory import SinglePatchRefinedSupportLagrangeTest as SinglePatchRefinedSupportLagrangeTest
from iga_test_factory import SinglePatchRefinedSupportNitscheTest as SinglePatchRefinedSupportNitscheTest
# Coupling/C_0 tests
from iga_test_factory import TwoPatchCouplingPenaltyShell3pTest as TwoPatchCouplingPenaltyShell3pTest
from iga_test_factory import TwoPatchCouplingLagrangeShell3pTest as TwoPatchCouplingLagrangeShell3pTest
from iga_test_factory import TwoPatchCouplingNitscheShell3pTest as TwoPatchCouplingNitscheShell3pTest
from iga_test_factory import TwoPatchRefinedCouplingPenaltyMembraneTest as TwoPatchRefinedCouplingPenaltyMembraneTest
from iga_test_factory import TwoPatchRefinedCouplingLagrangeMembraneTest as TwoPatchRefinedCouplingLagrangeMembraneTest
from iga_test_factory import TwoPatchRefinedCouplingNitscheMembraneTest as TwoPatchRefinedCouplingNitscheMembraneTest
# Rotation/G_1 coupling tests
from iga_test_factory import TwoPatchCantileverCouplingPenaltyTest as TwoPatchCantileverCouplingPenaltyTest
from iga_test_factory import TwoPatchCantileverRefinedCouplingPenaltyTest as TwoPatchCantileverRefinedCouplingPenaltyTest
# Nurbs Volume tests
from test_nurbs_volume_element import TestNurbsVolumeElement as TTestNurbsVolumeElements
# Modelers tests
from test_modelers import TestModelers as TTestModelers
from test_modelers_sbm import TestModelersSbm as TTestModelersSbm
from test_import_nurbs_modeler import TestImportNurbsModeler as TTestImportNurbsModeler
# Processes tests
from test_map_nurbs_volume_results_to_embedded_geometry_process import TestMapNurbsVolumeResultsToEmbeddedGeometryProcess as TTestMapNurbsVolumeResultsToEmbeddedGeometryProcess
# Fluid Element and Conditions tests
from applications.IgaApplication.tests.test_stokes_elements_and_conditions import FluidTests as TTestFluid
from applications.IgaApplication.tests.test_stokes_sbm_conditions import SbmStokesTests as TTestSbmStokes

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
        # Single patch test - checks iga essentials
        SinglePatchTest,
        # Truss tests
        TTrussElementTests,
        # Structural Elements tests
        TSolidIGAElementTests,
        # Sbm Elements tests
        TSbmSolidTests,
        # Membrane tests
        MembraneSinglePatchFourPointSailLinearStatic,
        MembraneSinglePatchFourPointSailNonLinearStatic,
        # 3p Shell KL
        TShell3pElementTests,
        ScordelisRoofShell3pTest,
        LinearBeamShell3pTest,
        LinearBeamShell3pAdditiveSchwarzTest,
        # 5p Shell Director
        ScordelisRoofShell5pTest,
        # Weak support tests
        SinglePatchRefinedSupportPenaltyTest,
        SinglePatchRefinedSupportLagrangeTest,
        # Coupling tests
        TwoPatchCouplingPenaltyShell3pTest,
        TwoPatchCouplingLagrangeShell3pTest,
        TwoPatchRefinedCouplingPenaltyMembraneTest,
        TwoPatchRefinedCouplingLagrangeMembraneTest,
        # Rotation/G_1 coupling tests
        TwoPatchCantileverCouplingPenaltyTest,
        TwoPatchCantileverRefinedCouplingPenaltyTest,
        # Volumes
        TTestNurbsVolumeElements,
        # Modelers
        TTestModelers,
        TTestModelersSbm,
        TTestMapNurbsVolumeResultsToEmbeddedGeometryProcess,
        # Fluids
        TTestFluid,
        TTestSbmStokes
    ]))

    if has_linear_solvers_application:
        from KratosMultiphysics import LinearSolversApplication
        if LinearSolversApplication.HasFEAST():
            smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([
                # Weak support Nitsche test
                SinglePatchRefinedSupportNitscheTest,
                # Coupling Nitsche tests
                TwoPatchCouplingNitscheShell3pTest,
                TwoPatchRefinedCouplingNitscheMembraneTest
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
        ]))

    validationSuite = suites['validation']
    validationSuite.addTests(nightSuite)
    validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([
        TShell5pHierarchicLinearScordelisTest,
        TShell5pHierarchicNonLinearThickBeamTest
    ]))

    allSuite = suites['all']
    allSuite.addTests(validationSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
