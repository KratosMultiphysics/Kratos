# import Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import subprocess
import subprocess

# Using kratos_utilities
import KratosMultiphysics.kratos_utilities as kratos_utilities
if kratos_utilities.CheckIfApplicationsAvailable("ExternalSolversApplication"):
    has_external_solvers_application = True
else:
    has_external_solvers_application = False

import os
import sys
def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)
sys.path.insert(0, GetFilePath('../../../kratos/tests/'))

# Import the tests o test_classes to create the suits
## SMALL TESTS
# Exact integration tests
from test_process_factory import TestProcessFactory as TTestProcessFactory
from test_double_curvature_integration import TestDoubleCurvatureIntegration as TTestDoubleCurvatureIntegration
from test_dynamic_search import TestDynamicSearch as TTestDynamicSearch
from test_mortar_mapper import TestMortarMapperCore as TTestMortarMapperCore

# Mesh tying tests
from SmallTests import SimplePatchTestTwoDMeshTying            as TSimplePatchTestTwoDMeshTying
from SmallTests import SimpleSlopePatchTestTwoDMeshTying       as TSimpleSlopePatchTestTwoDMeshTying
from SmallTests import SimplestPatchTestThreeDMeshTying        as TSimplestPatchTestThreeDMeshTying
from SmallTests import SimplestPatchTestThreeDTriQuadMeshTying as TSimplestPatchTestThreeDTriQuadMeshTying
from SmallTests import SimplestPatchTestThreeDQuadTriMeshTying as TSimplestPatchTestThreeDQuadTriMeshTying
from SmallTests import SimplePatchTestThreeDMeshTying          as TSimplePatchTestThreeDMeshTying

# ALM frictionless tests
from SmallTests import ALMHyperSimplePatchTestContact                              as TALMHyperSimplePatchTestContact
from SmallTests import ALMHyperSimplePatchTestWithEliminationContact               as TALMHyperSimplePatchTestWithEliminationContact
from SmallTests import ALMHyperSimplePatchTestWithEliminationWithConstraintContact as TALMHyperSimplePatchTestWithEliminationWithConstraintContact
from SmallTests import ALMHyperSimpleSlopePatchTestContact                         as TALMHyperSimpleSlopePatchTestContact
from SmallTests import ALMTwoDPatchComplexGeomTestContact                          as TALMTwoDPatchComplexGeomTestContact
from SmallTests import ALMTwoDPatchComplexGeomSlopeTestContact                     as TALMTwoDPatchComplexGeomSlopeTestContact
from SmallTests import ALMSimplePatchTestContact                                   as TALMSimplePatchTestContact
from SmallTests import ALMSimpleSlopePatchTestContact                              as TALMSimpleSlopePatchTestContact
from SmallTests import ALMSimplePatchNotMatchingATestContact                       as TALMSimplePatchNotMatchingATestContact
from SmallTests import ALMSimplePatchNotMatchingBTestContact                       as TALMSimplePatchNotMatchingBTestContact
from SmallTests import ALMThreeDSimplestPatchMatchingTestContact                   as TALMThreeDSimplestPatchMatchingTestContact
from SmallTests import ALMThreeDSimplestPatchTestTriQuadContact                    as TALMThreeDSimplestPatchTestTriQuadContact
from SmallTests import ALMThreeDSimplestPatchTestQuadTriContact                    as TALMThreeDSimplestPatchTestQuadTriContact
from SmallTests import ALMThreeDSimplestPatchMatchingAdaptativeTestContact         as TALMThreeDSimplestPatchMatchingAdaptativeTestContact
from SmallTests import ALMThreeDSimplestPatchMatchingSlopeTestContact              as TALMThreeDSimplestPatchMatchingSlopeTestContact
from SmallTests import ALMThreeDPatchComplexGeomTestContact                        as TALMThreeDPatchComplexGeomTestContact
from SmallTests import ALMThreeDPatchMatchingTestContact                           as TALMTThreeDPatchMatchingTestContact
from SmallTests import ALMThreeDPatchNotMatchingTestContact                        as TALMThreeDPatchNotMatchingTestContact

# Penalty frictionless tests
from SmallTests import PenaltyFrictionlessHyperSimplePatchFrictionalTestContact as TPenaltyFrictionlessHyperSimplePatchFrictionalTestContact
from SmallTests import PenaltyThreeDSimplestPatchMatchingTestContact as TPenaltyThreeDSimplestPatchMatchingTestContact
from NightlyTests import ExplicitPenaltyThreeDSimplestPatchMatchingTestContact as TExplicitPenaltyThreeDSimplestPatchMatchingTestContact

# Components ALM frictionless tests
from SmallTests import ComponentsALMHyperSimpleTrianglePatchTestContact                      as TComponentsALMHyperSimpleTrianglePatchTestContact
from SmallTests import ComponentsALMHyperSimplePatchTestContact                              as TComponentsALMHyperSimplePatchTestContact
from SmallTests import ComponentsALMHyperSimplePatchTestWithEliminationContact               as TComponentsALMHyperSimplePatchTestWithEliminationContact
from SmallTests import ComponentsALMHyperSimplePatchTestWithEliminationWithConstraintContact as TComponentsALMHyperSimplePatchTestWithEliminationWithConstraintContact
from SmallTests import ComponentsALMHyperSimpleSlopePatchTestContact                         as TComponentsALMHyperSimpleSlopePatchTestContact
from SmallTests import ComponentsALMTwoDPatchComplexGeomTestContact                          as TComponentsALMTwoDPatchComplexGeomTestContact
from SmallTests import ComponentsALMTwoDPatchComplexGeomSlopeTestContact                     as TComponentsALMTwoDPatchComplexGeomSlopeTestContact
from SmallTests import ComponentsALMSimplePatchTestContact                                   as TComponentsALMSimplePatchTestContact
from SmallTests import ComponentsALMSimpleSlopePatchTestContact                              as TComponentsALMSimpleSlopePatchTestContact
from SmallTests import ComponentsALMSimplePatchNotMatchingATestContact                       as TComponentsALMSimplePatchNotMatchingATestContact
from SmallTests import ComponentsALMSimplePatchNotMatchingBTestContact                       as TComponentsALMSimplePatchNotMatchingBTestContact
from SmallTests import ComponentsALMThreeDSimplestPatchMatchingTestContact                   as TComponentsALMThreeDSimplestPatchMatchingTestContact
from SmallTests import ComponentsALMThreeDSimplestPatchMatchingAdaptativeTestContact         as TComponentsALMThreeDSimplestPatchMatchingAdaptativeTestContact
from SmallTests import ComponentsALMThreeDSimplestPatchMatchingSlopeTestContact              as TComponentsALMThreeDSimplestPatchMatchingSlopeTestContact
from SmallTests import ComponentsALMThreeDPatchComplexGeomTestContact                        as TComponentsALMThreeDPatchComplexGeomTestContact
from SmallTests import ComponentsALMThreeDPatchMatchingTestContact                           as TComponentsALMTThreeDPatchMatchingTestContact
from SmallTests import ComponentsALMThreeDPatchNotMatchingTestContact                        as TComponentsALMThreeDPatchNotMatchingTestContact

# ALM frictional tests
from SmallTests import ALMHyperSimplePatchFrictionalTestContact                      as TALMHyperSimplePatchFrictionalTestContact

# Penalty frictional tests
from SmallTests import PenaltyNoFrictionHyperSimplePatchFrictionalTestContact        as TPenaltyNoFrictionHyperSimplePatchFrictionalTestContact
from SmallTests import PenaltyPerfectStickHyperSimplePatchFrictionalTestContact      as TPenaltyPerfectStickHyperSimplePatchFrictionalTestContact
from SmallTests import PenaltyThresholdSlipHyperSimplePatchFrictionalTestContact     as TPenaltyThresholdSlipHyperSimplePatchFrictionalTestContact
from SmallTests import PenaltyHyperSimplePatchFrictionalSlipTestContact              as TPenaltyHyperSimplePatchFrictionalSlipTestContact
from SmallTests import PenaltyHyperSimplePatchFrictionalStickTestContact             as TPenaltyHyperSimplePatchFrictionalStickTestContact

## NIGTHLY TESTS
# ALM frictionless tests
from NightlyTests import ALMTaylorPatchTestContact           as TALMTaylorPatchTestContact
from NightlyTests import ALMHertzSimpleTestContact           as TALMHertzSimpleTestContact
from NightlyTests import ALMHertzSimpleSphereTestContact     as TALMHertzSimpleSphereTestContact
from NightlyTests import ALMHertzSphereTestContact           as TALMHertzSphereTestContact
from NightlyTests import ALMHertzCompleteTestContact         as TALMHertzCompleteTestContact

# Components ALM frictionless tests
from NightlyTests import ComponentsALMTaylorPatchTestContact           as TComponentsALMTaylorPatchTestContact
from NightlyTests import ComponentsALMHertzSimpleTestContact           as TComponentsALMHertzSimpleTestContact
from NightlyTests import ComponentsALMHertzSimpleSphereTestContact     as TComponentsALMHertzSimpleSphereTestContact
from NightlyTests import ComponentsALMHertzSphereTestContact           as TComponentsALMHertzSphereTestContact
from NightlyTests import ComponentsALMHertzCompleteTestContact         as TComponentsALMHertzCompleteTestContact

# ALM frictionless tests
from NightlyTests import ALMTaylorPatchFrictionalTestContact           as TALMTaylorPatchFrictionalTestContact
from NightlyTests import ALMPureFrictionalTestContact                  as TALMPureFrictionalTestContact

## VALIDATION TESTS
from ValidationTests import LargeDisplacementPatchTestHexa as TLargeDisplacementPatchTestHexa

# ALM frictionless tests
from ValidationTests import ALMTaylorPatchDynamicTestContact as TALMTaylorPatchDynamicTestContact
from ValidationTests import ALMMeshMovingMatchingTestContact    as TALMMeshMovingMatchingTestContact
from ValidationTests import ALMMeshMovingNotMatchingTestContact as TALMMeshMovingNotMatchingTestContact
from ValidationTests import ALMIroningTestContact    as TALMIroningTestContact
from ValidationTests import ALMIroningDieTestContact as TALMIroningDieTestContact
from ValidationTests import ALMLargeDisplacementPatchTestTetra as TALMLargeDisplacementPatchTestTetra
from ValidationTests import ALMLargeDisplacementPatchTestHexa as TALMLargeDisplacementPatchTestHexa
from ValidationTests import ALMMultiLayerContactTest as TALMMultiLayerContactTest

# Components ALM frictionless tests
from ValidationTests import ComponentsALMTaylorPatchDynamicTestContact as TComponentsALMTaylorPatchDynamicTestContact
from ValidationTests import ComponentsALMMeshMovingMatchingTestContact    as TComponentsALMMeshMovingMatchingTestContact
from ValidationTests import ComponentsALMMeshMovingNotMatchingTestContact as TComponentsALMMeshMovingNotMatchingTestContact
from ValidationTests import ComponentsALMLargeDisplacementPatchTestTetra as TComponentsALMLargeDisplacementPatchTestTetra
from ValidationTests import ComponentsALMLargeDisplacementPatchTestHexa as TComponentsALMLargeDisplacementPatchTestHexa
from ValidationTests import ComponentsALMMultiLayerContactTest as TComponentsALMMultiLayerContactTest

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

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']

    # Test ProcessFactoryUtility
    smallSuite.addTest(TTestProcessFactory('test_process_factory'))
    smallSuite.addTest(TTestProcessFactory('test_processes_list_factory'))

    # Mesh tying tests
    smallSuite.addTest(TSimplePatchTestTwoDMeshTying('test_execution'))
    smallSuite.addTest(TSimpleSlopePatchTestTwoDMeshTying('test_execution'))
    smallSuite.addTest(TSimplestPatchTestThreeDMeshTying('test_execution'))
    smallSuite.addTest(TSimplestPatchTestThreeDTriQuadMeshTying('test_execution'))
    smallSuite.addTest(TSimplestPatchTestThreeDQuadTriMeshTying('test_execution'))
    smallSuite.addTest(TSimplePatchTestThreeDMeshTying('test_execution'))

    # ALM frictionless tests
    smallSuite.addTest(TALMHyperSimplePatchTestContact('test_execution'))
    smallSuite.addTest(TALMHyperSimplePatchTestWithEliminationContact('test_execution'))
    smallSuite.addTest(TALMHyperSimplePatchTestWithEliminationWithConstraintContact('test_execution'))
    smallSuite.addTest(TALMHyperSimpleSlopePatchTestContact('test_execution'))
    smallSuite.addTest(TALMTwoDPatchComplexGeomTestContact('test_execution'))
    smallSuite.addTest(TALMTwoDPatchComplexGeomSlopeTestContact('test_execution'))
    smallSuite.addTest(TALMSimplePatchTestContact('test_execution'))
    smallSuite.addTest(TALMSimpleSlopePatchTestContact('test_execution'))
    smallSuite.addTest(TALMSimplePatchNotMatchingATestContact('test_execution'))
    smallSuite.addTest(TALMSimplePatchNotMatchingBTestContact('test_execution'))
    smallSuite.addTest(TALMThreeDSimplestPatchMatchingTestContact('test_execution'))
    smallSuite.addTest(TALMThreeDSimplestPatchTestTriQuadContact('test_execution'))
    smallSuite.addTest(TALMThreeDSimplestPatchTestQuadTriContact('test_execution'))
    smallSuite.addTest(TALMThreeDSimplestPatchMatchingAdaptativeTestContact('test_execution'))
    smallSuite.addTest(TALMThreeDSimplestPatchMatchingSlopeTestContact('test_execution'))
    smallSuite.addTest(TALMThreeDPatchComplexGeomTestContact('test_execution'))

    # Penalty frictionless tests
    smallSuite.addTest(TPenaltyFrictionlessHyperSimplePatchFrictionalTestContact('test_execution'))
    smallSuite.addTest(TPenaltyThreeDSimplestPatchMatchingTestContact('test_execution'))

    # Components ALM frictionless tests
    smallSuite.addTest(TComponentsALMHyperSimpleTrianglePatchTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMHyperSimplePatchTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMHyperSimplePatchTestWithEliminationContact('test_execution'))
    smallSuite.addTest(TComponentsALMHyperSimplePatchTestWithEliminationWithConstraintContact('test_execution'))
    smallSuite.addTest(TComponentsALMHyperSimpleSlopePatchTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMTwoDPatchComplexGeomTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMTwoDPatchComplexGeomSlopeTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMSimplePatchTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMSimpleSlopePatchTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMSimplePatchNotMatchingATestContact('test_execution'))
    smallSuite.addTest(TComponentsALMSimplePatchNotMatchingBTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMThreeDSimplestPatchMatchingTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMThreeDSimplestPatchMatchingAdaptativeTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMThreeDSimplestPatchMatchingSlopeTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMThreeDPatchComplexGeomTestContact('test_execution'))

    # ALM frictional tests
    smallSuite.addTest(TALMHyperSimplePatchFrictionalTestContact('test_execution'))

    # Penalty frictional tests
    smallSuite.addTest(TPenaltyNoFrictionHyperSimplePatchFrictionalTestContact('test_execution'))
    smallSuite.addTest(TPenaltyPerfectStickHyperSimplePatchFrictionalTestContact('test_execution'))
    smallSuite.addTest(TPenaltyThresholdSlipHyperSimplePatchFrictionalTestContact('test_execution'))
    smallSuite.addTest(TPenaltyHyperSimplePatchFrictionalSlipTestContact('test_execution'))
    smallSuite.addTest(TPenaltyHyperSimplePatchFrictionalStickTestContact('test_execution'))

    # Fill with all small tests
    nightSuite.addTests(smallSuite)

    # Exact integration tests
    nightSuite.addTest(TTestDoubleCurvatureIntegration('test_moving_mesh_integration_quad'))

    # Mortar mapping
    nightSuite.addTest(TTestMortarMapperCore('test_less_basic_mortar_mapping_triangle'))
    nightSuite.addTest(TTestMortarMapperCore('test_simple_curvature_mortar_mapping_triangle'))

    # ALM frictionless tests
    nightSuite.addTest(TALMTThreeDPatchMatchingTestContact('test_execution'))
    nightSuite.addTest(TALMThreeDPatchNotMatchingTestContact('test_execution'))
    nightSuite.addTest(TALMTaylorPatchTestContact('test_execution'))
    nightSuite.addTest(TALMHertzSimpleSphereTestContact('test_execution'))

    # Components ALM frictionless tests
    nightSuite.addTest(TComponentsALMTThreeDPatchMatchingTestContact('test_execution'))
    nightSuite.addTest(TComponentsALMThreeDPatchNotMatchingTestContact('test_execution'))
    nightSuite.addTest(TComponentsALMTaylorPatchTestContact('test_execution'))
    nightSuite.addTest(TComponentsALMHertzSimpleSphereTestContact('test_execution'))

    # ALM frictional tests
    nightSuite.addTest(TALMTaylorPatchFrictionalTestContact('test_execution'))
    nightSuite.addTest(TALMPureFrictionalTestContact('test_execution'))

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(nightSuite)

    # ALM frictionless tests
    #nightSuite.addTest(TALMHertzSphereTestContact('test_execution'))
    validationSuite.addTest(TALMHertzSimpleTestContact('test_execution'))
    validationSuite.addTest(TALMHertzCompleteTestContact('test_execution'))

    # Components ALM frictionless tests
    #nightSuite.addTest(TComponentsALMHertzSphereTestContact('test_execution'))
    validationSuite.addTest(TComponentsALMHertzSimpleTestContact('test_execution'))
    validationSuite.addTest(TComponentsALMHertzCompleteTestContact('test_execution'))

    # Penalty frictionless tests
    validationSuite.addTest(TExplicitPenaltyThreeDSimplestPatchMatchingTestContact('test_execution'))

    # Exact integration tests
    validationSuite.addTest(TTestDoubleCurvatureIntegration('test_double_curvature_integration_triangle'))
    validationSuite.addTest(TTestDoubleCurvatureIntegration('test_double_curvature_integration_quad'))
    validationSuite.addTest(TTestDoubleCurvatureIntegration('test_moving_mesh_integration_quad'))

    # Dynamic search
    validationSuite.addTest(TTestDynamicSearch('test_dynamic_search_triangle'))
    validationSuite.addTest(TTestDynamicSearch('test_dynamic_search_quad'))

    # Mortar mapping
    validationSuite.addTest(TTestMortarMapperCore('test_mortar_mapping_triangle'))
    validationSuite.addTest(TTestMortarMapperCore('test_mortar_mapping_quad'))

    validationSuite.addTest(TLargeDisplacementPatchTestHexa('test_execution'))

    # ALM frictionless tests
    validationSuite.addTest(TALMTaylorPatchDynamicTestContact('test_execution'))
    validationSuite.addTest(TALMMeshMovingMatchingTestContact('test_execution'))
    validationSuite.addTest(TALMMeshMovingNotMatchingTestContact('test_execution'))
    #validationSuite.addTest(TALMIroningTestContact('test_execution'))
    #validationSuite.addTest(TALMIroningDieTestContact('test_execution'))
    validationSuite.addTest(TALMLargeDisplacementPatchTestTetra('test_execution'))
    validationSuite.addTest(TALMLargeDisplacementPatchTestHexa('test_execution'))
    validationSuite.addTest(TALMMultiLayerContactTest('test_execution'))

    # Components ALM frictionless tests
    validationSuite.addTest(TComponentsALMTaylorPatchDynamicTestContact('test_execution'))
    validationSuite.addTest(TComponentsALMMeshMovingMatchingTestContact('test_execution'))
    validationSuite.addTest(TComponentsALMMeshMovingNotMatchingTestContact('test_execution'))
    validationSuite.addTest(TComponentsALMLargeDisplacementPatchTestTetra('test_execution'))
    validationSuite.addTest(TComponentsALMLargeDisplacementPatchTestHexa('test_execution'))
    validationSuite.addTest(TComponentsALMMultiLayerContactTest('test_execution'))

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            ### STANDALONE
            TTestProcessFactory,
            TTestDoubleCurvatureIntegration,
            TTestDynamicSearch,
            TTestMortarMapperCore,
            ### SMALL
            TSimplePatchTestTwoDMeshTying,
            TSimpleSlopePatchTestTwoDMeshTying,
            TSimplestPatchTestThreeDMeshTying,
            TSimplestPatchTestThreeDTriQuadMeshTying,
            TSimplestPatchTestThreeDQuadTriMeshTying,
            TSimplePatchTestThreeDMeshTying,
            TALMHyperSimplePatchTestContact,
            TALMHyperSimplePatchTestWithEliminationContact,
            TALMHyperSimplePatchTestWithEliminationWithConstraintContact,
            TALMHyperSimpleSlopePatchTestContact,
            TALMTwoDPatchComplexGeomTestContact,
            TALMTwoDPatchComplexGeomSlopeTestContact,
            TALMSimplePatchTestContact,
            TALMSimpleSlopePatchTestContact,
            TALMSimplePatchNotMatchingATestContact,
            TALMSimplePatchNotMatchingBTestContact,
            TALMThreeDSimplestPatchMatchingTestContact,
            TALMThreeDSimplestPatchTestTriQuadContact,
            TALMThreeDSimplestPatchTestQuadTriContact,
            TALMThreeDSimplestPatchMatchingAdaptativeTestContact,
            TALMThreeDSimplestPatchMatchingSlopeTestContact,
            TALMThreeDPatchComplexGeomTestContact,
            TALMTThreeDPatchMatchingTestContact,
            TALMThreeDPatchNotMatchingTestContact,
            TPenaltyFrictionlessHyperSimplePatchFrictionalTestContact,
            TPenaltyThreeDSimplestPatchMatchingTestContact,
            TExplicitPenaltyThreeDSimplestPatchMatchingTestContact,
            TComponentsALMHyperSimpleTrianglePatchTestContact,
            TComponentsALMHyperSimplePatchTestContact,
            TComponentsALMHyperSimplePatchTestWithEliminationContact,
            TComponentsALMHyperSimplePatchTestWithEliminationWithConstraintContact,
            TComponentsALMHyperSimpleSlopePatchTestContact,
            TComponentsALMTwoDPatchComplexGeomTestContact,
            TComponentsALMTwoDPatchComplexGeomSlopeTestContact,
            TComponentsALMSimplePatchTestContact,
            TComponentsALMSimpleSlopePatchTestContact,
            TComponentsALMSimplePatchNotMatchingATestContact,
            TComponentsALMSimplePatchNotMatchingBTestContact,
            TComponentsALMThreeDSimplestPatchMatchingTestContact,
            TComponentsALMThreeDSimplestPatchMatchingAdaptativeTestContact,
            TComponentsALMThreeDSimplestPatchMatchingSlopeTestContact,
            TComponentsALMThreeDPatchComplexGeomTestContact,
            TComponentsALMTThreeDPatchMatchingTestContact,
            TComponentsALMThreeDPatchNotMatchingTestContact,
            TALMHyperSimplePatchFrictionalTestContact,
            TPenaltyNoFrictionHyperSimplePatchFrictionalTestContact,
            TPenaltyPerfectStickHyperSimplePatchFrictionalTestContact,
            TPenaltyThresholdSlipHyperSimplePatchFrictionalTestContact,
            TPenaltyHyperSimplePatchFrictionalSlipTestContact,
            TPenaltyHyperSimplePatchFrictionalStickTestContact,
            #### NIGTHLY
            TALMTaylorPatchTestContact,
            TALMHertzSimpleTestContact,
            TALMHertzSimpleSphereTestContact,
            #####TALMHertzSphereTestContact,  # FIXME: This test requieres the axisymmetric to work (memmory error, correct it)
            TALMHertzCompleteTestContact,
            TComponentsALMTaylorPatchTestContact,
            TComponentsALMHertzSimpleTestContact,
            TComponentsALMHertzSimpleSphereTestContact,
            #####TComponentsALMHertzSphereTestContact,  # FIXME: This test requieres the axisymmetric to work (memmory error, correct it)
            TComponentsALMHertzCompleteTestContact,
            TALMTaylorPatchFrictionalTestContact,
            TALMPureFrictionalTestContact,
            #### VALIDATION
            TALMTaylorPatchDynamicTestContact,
            TALMMeshMovingMatchingTestContact,
            TALMMeshMovingNotMatchingTestContact,
            ##TALMIroningTestContact,
            ##TALMIroningDieTestContact,
            TLargeDisplacementPatchTestHexa,
            TALMLargeDisplacementPatchTestTetra,
            TALMLargeDisplacementPatchTestHexa,
            TALMMultiLayerContactTest,
            TComponentsALMTaylorPatchDynamicTestContact,
            TComponentsALMMeshMovingMatchingTestContact,
            TComponentsALMMeshMovingNotMatchingTestContact,
            TComponentsALMLargeDisplacementPatchTestTetra,
            TComponentsALMLargeDisplacementPatchTestHexa,
            TComponentsALMMultiLayerContactTest,
        ])
    )

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
