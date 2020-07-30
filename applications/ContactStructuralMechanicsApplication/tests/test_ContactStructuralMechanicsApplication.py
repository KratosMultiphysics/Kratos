# import Kratos
import KratosMultiphysics
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

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
from test_process_factory import TestProcessFactory                          as TTestProcessFactory
from test_check_normals_process import TestCheckNormals                      as TTestCheckNormals
from test_double_curvature_integration import TestDoubleCurvatureIntegration as TTestDoubleCurvatureIntegration
from test_dynamic_search import TestDynamicSearch                            as TTestDynamicSearch
from test_mortar_mapper import TestMortarMapperCore                          as TTestMortarMapperCore

# Mesh tying tests
from SmallTests import SimplePatchTestTwoDMeshTying            as TSimplePatchTestTwoDMeshTying
from SmallTests import SimpleSlopePatchTestTwoDMeshTying       as TSimpleSlopePatchTestTwoDMeshTying
from SmallTests import SimplestPatchTestThreeDMeshTying        as TSimplestPatchTestThreeDMeshTying

# ALM frictionless tests
from SmallTests import ALMHyperSimplePatchTestContact                              as TALMHyperSimplePatchTestContact
from SmallTests import ALMHyperSimplePatchTrianglesTestContact                     as TALMHyperSimplePatchTrianglesTestContact
from SmallTests import ALMHyperSimplePatchTestWithEliminationContact               as TALMHyperSimplePatchTestWithEliminationContact
from SmallTests import ALMHyperSimplePatchTestWithEliminationWithConstraintContact as TALMHyperSimplePatchTestWithEliminationWithConstraintContact
from SmallTests import ALMHyperSimpleSlopePatchTestContact                         as TALMHyperSimpleSlopePatchTestContact
from SmallTests import ALMThreeDSimplestPatchMatchingTestContact                   as TALMThreeDSimplestPatchMatchingTestContact

# Penalty frictionless tests
from SmallTests import PenaltyFrictionlessHyperSimplePatchFrictionalTestContact as TPenaltyFrictionlessHyperSimplePatchFrictionalTestContact
from SmallTests import PenaltyThreeDSimplestPatchMatchingTestContact            as TPenaltyThreeDSimplestPatchMatchingTestContact

# Components ALM frictionless tests
from SmallTests import ComponentsALMHyperSimpleTrianglePatchTestContact                      as TComponentsALMHyperSimpleTrianglePatchTestContact
from SmallTests import ComponentsALMHyperSimplePatchTestContact                              as TComponentsALMHyperSimplePatchTestContact
from SmallTests import ComponentsALMHyperSimplePatchTestWithEliminationContact               as TComponentsALMHyperSimplePatchTestWithEliminationContact
from SmallTests import ComponentsALMHyperSimplePatchTestWithEliminationWithConstraintContact as TComponentsALMHyperSimplePatchTestWithEliminationWithConstraintContact
from SmallTests import ComponentsALMHyperSimpleSlopePatchTestContact                         as TComponentsALMHyperSimpleSlopePatchTestContact
from SmallTests import ComponentsALMThreeDSimplestPatchMatchingTestContact                   as TComponentsALMThreeDSimplestPatchMatchingTestContact

# ALM frictional tests
from SmallTests import ALMHyperSimplePatchFrictionalTestContact                      as TALMHyperSimplePatchFrictionalTestContact
from SmallTests import ALMNoFrictionHyperSimplePatchFrictionalTestContact            as TALMNoFrictionHyperSimplePatchFrictionalTestContact
from SmallTests import ALMPerfectStickHyperSimplePatchFrictionalTestContact          as TALMPerfectStickHyperSimplePatchFrictionalTestContact
from SmallTests import ALMThresholdSlipHyperSimplePatchFrictionalTestContact         as TALMThresholdSlipHyperSimplePatchFrictionalTestContact
from SmallTests import ALMHyperSimplePatchFrictionalSlipTestContact                  as TALMHyperSimplePatchFrictionalSlipTestContact
from SmallTests import ALMHyperSimplePatchFrictionalStickTestContact                 as TALMHyperSimplePatchFrictionalStickTestContact

# Penalty frictional tests
from SmallTests import PenaltyNoFrictionHyperSimplePatchFrictionalTestContact        as TPenaltyNoFrictionHyperSimplePatchFrictionalTestContact
from SmallTests import PenaltyPerfectStickHyperSimplePatchFrictionalTestContact      as TPenaltyPerfectStickHyperSimplePatchFrictionalTestContact
from SmallTests import PenaltyThresholdSlipHyperSimplePatchFrictionalTestContact     as TPenaltyThresholdSlipHyperSimplePatchFrictionalTestContact
from SmallTests import PenaltyHyperSimplePatchFrictionalSlipTestContact              as TPenaltyHyperSimplePatchFrictionalSlipTestContact
from SmallTests import PenaltyHyperSimplePatchFrictionalStickTestContact             as TPenaltyHyperSimplePatchFrictionalStickTestContact

# MPC Contact tests
from SmallTests import TwoDSimplestPatchMatchingTestContact               as TTwoDSimplestPatchMatchingTestContact
from SmallTests import TwoDSimplestWithFrictionPatchMatchingTestContact   as TTwoDSimplestWithFrictionPatchMatchingTestContact
from SmallTests import ThreeDSimplestPatchMatchingTestContact             as TThreeDSimplestPatchMatchingTestContact
from SmallTests import ThreeDSimplestWithFrictionPatchMatchingTestContact as TThreeDSimplestWithFrictionPatchMatchingTestContact

## NIGTHLY TESTS
# Mesh tying tests
from NightlyTests import SimplestPatchTestThreeDTriQuadMeshTying as TSimplestPatchTestThreeDTriQuadMeshTying
from NightlyTests import SimplestPatchTestThreeDQuadTriMeshTying as TSimplestPatchTestThreeDQuadTriMeshTying
from NightlyTests import SimplePatchTestThreeDMeshTying          as TSimplePatchTestThreeDMeshTying

# ALM frictionless tests
from NightlyTests import ALMTwoDPatchComplexGeomTestContact                          as TALMTwoDPatchComplexGeomTestContact
from NightlyTests import ALMTwoDPatchComplexGeomSlopeTestContact                     as TALMTwoDPatchComplexGeomSlopeTestContact
from NightlyTests import ALMSimplePatchTestContact                                   as TALMSimplePatchTestContact
from NightlyTests import ALMSimpleSlopePatchTestContact                              as TALMSimpleSlopePatchTestContact
from NightlyTests import ALMSimplePatchNotMatchingATestContact                       as TALMSimplePatchNotMatchingATestContact
from NightlyTests import ALMSimplePatchNotMatchingBTestContact                       as TALMSimplePatchNotMatchingBTestContact
from NightlyTests import ALMThreeDSimplestPatchTestTriQuadContact                    as TALMThreeDSimplestPatchTestTriQuadContact
from NightlyTests import ALMThreeDSimplestPatchTestQuadTriContact                    as TALMThreeDSimplestPatchTestQuadTriContact
from NightlyTests import ALMThreeDSimplestPatchMatchingAdaptativeTestContact         as TALMThreeDSimplestPatchMatchingAdaptativeTestContact
from NightlyTests import ALMThreeDSimplestPatchMatchingSlopeTestContact              as TALMThreeDSimplestPatchMatchingSlopeTestContact
from NightlyTests import ALMThreeDPatchComplexGeomTestContact                        as TALMThreeDPatchComplexGeomTestContact
from NightlyTests import ALMThreeDPatchMatchingTestContact                           as TALMTThreeDPatchMatchingTestContact
from NightlyTests import ALMThreeDPatchNotMatchingTestContact                        as TALMThreeDPatchNotMatchingTestContact
from NightlyTests import ALMTaylorPatchTestContact                                   as TALMTaylorPatchTestContact
from NightlyTests import ALMHertzSimpleTestContact                                   as TALMHertzSimpleTestContact
from NightlyTests import ALMHertzSimpleSphereTestContact                             as TALMHertzSimpleSphereTestContact
#from NightlyTests import ALMHertzSphereTestContact                                    as TALMHertzSphereTestContact
from NightlyTests import ALMHertzCompleteTestContact                                 as TALMHertzCompleteTestContact

# Components ALM frictionless tests
from NightlyTests import ComponentsALMTwoDPatchComplexGeomTestContact                          as TComponentsALMTwoDPatchComplexGeomTestContact
from NightlyTests import ComponentsALMTwoDPatchComplexGeomSlopeTestContact                     as TComponentsALMTwoDPatchComplexGeomSlopeTestContact
from NightlyTests import ComponentsALMSimplePatchTestContact                                   as TComponentsALMSimplePatchTestContact
from NightlyTests import ComponentsALMSimpleSlopePatchTestContact                              as TComponentsALMSimpleSlopePatchTestContact
from NightlyTests import ComponentsALMSimplePatchNotMatchingATestContact                       as TComponentsALMSimplePatchNotMatchingATestContact
from NightlyTests import ComponentsALMSimplePatchNotMatchingBTestContact                       as TComponentsALMSimplePatchNotMatchingBTestContact
from NightlyTests import ComponentsALMThreeDSimplestPatchTestTriQuadContact                    as TComponentsALMThreeDSimplestPatchTestTriQuadContact
from NightlyTests import ComponentsALMThreeDSimplestPatchTestQuadTriContact                    as TComponentsALMThreeDSimplestPatchTestQuadTriContact
from NightlyTests import ComponentsALMThreeDSimplestPatchMatchingAdaptativeTestContact         as TComponentsALMThreeDSimplestPatchMatchingAdaptativeTestContact
from NightlyTests import ComponentsALMThreeDSimplestPatchMatchingSlopeTestContact              as TComponentsALMThreeDSimplestPatchMatchingSlopeTestContact
from NightlyTests import ComponentsALMThreeDPatchComplexGeomTestContact                        as TComponentsALMThreeDPatchComplexGeomTestContact
from NightlyTests import ComponentsALMThreeDPatchMatchingTestContact                           as TComponentsALMTThreeDPatchMatchingTestContact
from NightlyTests import ComponentsALMThreeDPatchNotMatchingTestContact                        as TComponentsALMThreeDPatchNotMatchingTestContact
from NightlyTests import ComponentsALMTaylorPatchTestContact                                   as TComponentsALMTaylorPatchTestContact
from NightlyTests import ComponentsALMHertzSimpleTestContact                                   as TComponentsALMHertzSimpleTestContact
from NightlyTests import ComponentsALMHertzSimpleSphereTestContact                             as TComponentsALMHertzSimpleSphereTestContact
#from NightlyTests import ComponentsALMHertzSphereTestContact                                    as TComponentsALMHertzSphereTestContact
from NightlyTests import ComponentsALMHertzCompleteTestContact                                 as TComponentsALMHertzCompleteTestContact

# ALM frictional tests
from NightlyTests import ALMPureFrictionalTestContact                  as TALMPureFrictionalTestContact
from NightlyTests import ALMBasicFrictionTestContact                   as TALMBasicFrictionTestContact
from NightlyTests import ALMStaticEvolutionLoadFrictionTestContact     as TALMStaticEvolutionLoadFrictionTestContact
from NightlyTests import ALMEvolutionLoadFrictionTestContact           as TALMEvolutionLoadFrictionTestContact

# MPC Contact tests
from NightlyTests import ThreeDSimplestPatchMatchingSlopeTestContact as TThreeDSimplestPatchMatchingSlopeTestContact
from NightlyTests import ThreeDPatchMatchingTestContact              as TThreeDPatchMatchingTestContact
from NightlyTests import ThreeDPatchNotMatchingTestContact           as TThreeDPatchNotMatchingTestContact
from NightlyTests import BeamAxilSimpleContactTest                   as TBeamAxilSimpleContactTest
from NightlyTests import BeamAxilContactTest                         as TBeamAxilContactTest
from NightlyTests import BeamAxilTetraContactTest                    as TBeamAxilTetraContactTest
from NightlyTests import BeamContactTest                             as TBeamContactTest
from NightlyTests import BeamContactWithTyingTest                    as TBeamContactWithTyingTest
from NightlyTests import BeamContactWithFrictionTest                 as TBeamContactWithFrictionTest
from NightlyTests import PlateTest                                   as TPlateTest

## VALIDATION TESTS
from ValidationTests import LargeDisplacementPatchTestHexa as TLargeDisplacementPatchTestHexa
from ValidationTests import MeshTyingValidationTest        as TMeshTyingValidationTest

# ALM frictionless tests
from ValidationTests import ALMTaylorPatchDynamicTestContact    as TALMTaylorPatchDynamicTestContact
from ValidationTests import ALMMeshMovingMatchingTestContact    as TALMMeshMovingMatchingTestContact
from ValidationTests import ALMMeshMovingNotMatchingTestContact as TALMMeshMovingNotMatchingTestContact
#from ValidationTests import ALMIroningTestContact                as TALMIroningTestContact
#from ValidationTests import ALMIroningDieTestContact             as TALMIroningDieTestContact
from ValidationTests import ALMLargeDisplacementPatchTestTetra  as TALMLargeDisplacementPatchTestTetra
from ValidationTests import ALMLargeDisplacementPatchTestHexa   as TALMLargeDisplacementPatchTestHexa
from ValidationTests import ALMMultiLayerContactTest            as TALMMultiLayerContactTest
from ValidationTests import ALMSelfContactContactTest           as TALMSelfContactContactTest

# Penalty frictionless tests
from ValidationTests import ExplicitPenaltyThreeDSimplestPatchMatchingTestContact as TExplicitPenaltyThreeDSimplestPatchMatchingTestContact

# Components ALM frictionless tests
from ValidationTests import ComponentsALMTaylorPatchDynamicTestContact    as TComponentsALMTaylorPatchDynamicTestContact
from ValidationTests import ComponentsALMMeshMovingMatchingTestContact    as TComponentsALMMeshMovingMatchingTestContact
from ValidationTests import ComponentsALMMeshMovingNotMatchingTestContact as TComponentsALMMeshMovingNotMatchingTestContact
from ValidationTests import ComponentsALMLargeDisplacementPatchTestTetra  as TComponentsALMLargeDisplacementPatchTestTetra
from ValidationTests import ComponentsALMLargeDisplacementPatchTestHexa   as TComponentsALMLargeDisplacementPatchTestHexa
from ValidationTests import ComponentsALMMultiLayerContactTest            as TComponentsALMMultiLayerContactTest
from ValidationTests import ComponentsALMSelfContactContactTest           as TComponentsALMSelfContactContactTest

# ALM frictional tests
from ValidationTests import ALMTaylorPatchFrictionalTestContact                   as TALMTaylorPatchFrictionalTestContact
from ValidationTests import ALMMeshMovingMatchingTestFrictionalPureSlipContact    as TALMMeshMovingMatchingTestFrictionalPureSlipContact
from ValidationTests import ALMMeshMovingNotMatchingTestFrictionalPureSlipContact as TALMMeshMovingNotMatchingTestFrictionalPureSlipContact
from ValidationTests import ALMHertzTestFrictionalContact                         as TALMHertzTestFrictionalContact
from ValidationTests import ALMBlockTestFrictionalContact                         as TALMBlockTestFrictionalContact

##### VALIDATION TESTS #####
#from ValidationTests import MultiLayerContactTest as TMultiLayerContactTest

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
    nightlySuite = suites['nightly']

    ### BEGIN SMALL SUITE ###

    # Test ProcessFactoryUtility
    smallSuite.addTest(TTestProcessFactory('test_process_factory'))
    smallSuite.addTest(TTestCheckNormals('test_check_normals'))
    smallSuite.addTest(TTestCheckNormals('test_check_normals_quads'))
    smallSuite.addTest(TTestProcessFactory('test_processes_list_factory'))

    # Mesh tying tests
    smallSuite.addTest(TSimplePatchTestTwoDMeshTying('test_execution'))
    smallSuite.addTest(TSimpleSlopePatchTestTwoDMeshTying('test_execution'))
    smallSuite.addTest(TSimplestPatchTestThreeDMeshTying('test_execution'))

    # ALM frictionless tests
    smallSuite.addTest(TALMHyperSimplePatchTestContact('test_execution'))
    smallSuite.addTest(TALMHyperSimplePatchTrianglesTestContact('test_execution'))
    smallSuite.addTest(TALMHyperSimplePatchTestWithEliminationContact('test_execution'))
    smallSuite.addTest(TALMHyperSimplePatchTestWithEliminationWithConstraintContact('test_execution'))
    smallSuite.addTest(TALMHyperSimpleSlopePatchTestContact('test_execution'))
    smallSuite.addTest(TALMThreeDSimplestPatchMatchingTestContact('test_execution'))

    # Penalty frictionless tests
    smallSuite.addTest(TPenaltyFrictionlessHyperSimplePatchFrictionalTestContact('test_execution'))
    smallSuite.addTest(TPenaltyThreeDSimplestPatchMatchingTestContact('test_execution'))

    # Components ALM frictionless tests
    smallSuite.addTest(TComponentsALMHyperSimpleTrianglePatchTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMHyperSimplePatchTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMHyperSimplePatchTestWithEliminationContact('test_execution'))
    smallSuite.addTest(TComponentsALMHyperSimplePatchTestWithEliminationWithConstraintContact('test_execution'))
    smallSuite.addTest(TComponentsALMHyperSimpleSlopePatchTestContact('test_execution'))
    smallSuite.addTest(TComponentsALMThreeDSimplestPatchMatchingTestContact('test_execution'))

    # ALM frictional tests
    smallSuite.addTest(TALMHyperSimplePatchFrictionalTestContact('test_execution'))
    smallSuite.addTest(TALMNoFrictionHyperSimplePatchFrictionalTestContact('test_execution'))
    smallSuite.addTest(TALMPerfectStickHyperSimplePatchFrictionalTestContact('test_execution'))
    smallSuite.addTest(TALMThresholdSlipHyperSimplePatchFrictionalTestContact('test_execution'))
    smallSuite.addTest(TALMHyperSimplePatchFrictionalSlipTestContact('test_execution'))
    smallSuite.addTest(TALMHyperSimplePatchFrictionalStickTestContact('test_execution'))

    # Penalty frictional tests
    smallSuite.addTest(TPenaltyNoFrictionHyperSimplePatchFrictionalTestContact('test_execution'))
    smallSuite.addTest(TPenaltyPerfectStickHyperSimplePatchFrictionalTestContact('test_execution'))
    smallSuite.addTest(TPenaltyThresholdSlipHyperSimplePatchFrictionalTestContact('test_execution'))
    smallSuite.addTest(TPenaltyHyperSimplePatchFrictionalSlipTestContact('test_execution'))
    smallSuite.addTest(TPenaltyHyperSimplePatchFrictionalStickTestContact('test_execution'))

    # MPC contact test
    smallSuite.addTest(TTwoDSimplestPatchMatchingTestContact('test_execution'))
    smallSuite.addTest(TTwoDSimplestWithFrictionPatchMatchingTestContact('test_execution'))
    smallSuite.addTest(TThreeDSimplestPatchMatchingTestContact('test_execution'))
    smallSuite.addTest(TThreeDSimplestWithFrictionPatchMatchingTestContact('test_execution'))
    smallSuite.addTest(TThreeDSimplestPatchMatchingSlopeTestContact('test_execution'))
    smallSuite.addTest(TThreeDPatchMatchingTestContact('test_execution'))
    smallSuite.addTest(TThreeDPatchNotMatchingTestContact('test_execution'))

    ### END SMALL SUITE ###

    ### BEGIN NIGHTLY SUITE ###

    # Fill with all small tests
    nightlySuite.addTests(smallSuite)

    # Exact integration tests
    nightlySuite.addTest(TTestDoubleCurvatureIntegration('test_moving_mesh_integration_quad'))

    # Mortar mapping
    nightlySuite.addTest(TTestMortarMapperCore('test_less_basic_mortar_mapping_triangle'))
    nightlySuite.addTest(TTestMortarMapperCore('test_simple_curvature_mortar_mapping_triangle'))

    # Mesh tying tests
    nightlySuite.addTest(TSimplestPatchTestThreeDTriQuadMeshTying('test_execution'))
    nightlySuite.addTest(TSimplestPatchTestThreeDQuadTriMeshTying('test_execution'))
    nightlySuite.addTest(TSimplePatchTestThreeDMeshTying('test_execution'))

    # ALM frictionless tests
    nightlySuite.addTest(TALMTwoDPatchComplexGeomTestContact('test_execution'))
    nightlySuite.addTest(TALMTwoDPatchComplexGeomSlopeTestContact('test_execution'))
    nightlySuite.addTest(TALMSimplePatchTestContact('test_execution'))
    nightlySuite.addTest(TALMSimpleSlopePatchTestContact('test_execution'))
    nightlySuite.addTest(TALMSimplePatchNotMatchingATestContact('test_execution'))
    nightlySuite.addTest(TALMSimplePatchNotMatchingBTestContact('test_execution'))
    nightlySuite.addTest(TALMThreeDSimplestPatchTestTriQuadContact('test_execution'))
    nightlySuite.addTest(TALMThreeDSimplestPatchTestQuadTriContact('test_execution'))
    nightlySuite.addTest(TALMThreeDSimplestPatchMatchingAdaptativeTestContact('test_execution'))
    nightlySuite.addTest(TALMThreeDSimplestPatchMatchingSlopeTestContact('test_execution'))
    nightlySuite.addTest(TALMThreeDPatchComplexGeomTestContact('test_execution'))
    nightlySuite.addTest(TALMTThreeDPatchMatchingTestContact('test_execution'))
    nightlySuite.addTest(TALMThreeDPatchNotMatchingTestContact('test_execution'))
    nightlySuite.addTest(TALMTaylorPatchTestContact('test_execution'))
    nightlySuite.addTest(TALMHertzSimpleSphereTestContact('test_execution'))

    # Components ALM frictionless tests
    nightlySuite.addTest(TComponentsALMTwoDPatchComplexGeomTestContact('test_execution'))
    nightlySuite.addTest(TComponentsALMTwoDPatchComplexGeomSlopeTestContact('test_execution'))
    nightlySuite.addTest(TComponentsALMSimplePatchTestContact('test_execution'))
    nightlySuite.addTest(TComponentsALMSimpleSlopePatchTestContact('test_execution'))
    nightlySuite.addTest(TComponentsALMSimplePatchNotMatchingATestContact('test_execution'))
    nightlySuite.addTest(TComponentsALMSimplePatchNotMatchingBTestContact('test_execution'))
    nightlySuite.addTest(TComponentsALMThreeDSimplestPatchTestTriQuadContact('test_execution'))
    nightlySuite.addTest(TComponentsALMThreeDSimplestPatchTestQuadTriContact('test_execution'))
    nightlySuite.addTest(TComponentsALMThreeDSimplestPatchMatchingAdaptativeTestContact('test_execution'))
    nightlySuite.addTest(TComponentsALMThreeDSimplestPatchMatchingSlopeTestContact('test_execution'))
    nightlySuite.addTest(TComponentsALMThreeDPatchComplexGeomTestContact('test_execution'))
    nightlySuite.addTest(TComponentsALMTThreeDPatchMatchingTestContact('test_execution'))
    nightlySuite.addTest(TComponentsALMThreeDPatchNotMatchingTestContact('test_execution'))
    nightlySuite.addTest(TComponentsALMTaylorPatchTestContact('test_execution'))
    nightlySuite.addTest(TComponentsALMHertzSimpleSphereTestContact('test_execution'))

    # ALM frictional tests
    nightlySuite.addTest(TALMPureFrictionalTestContact('test_execution'))
    nightlySuite.addTest(TALMBasicFrictionTestContact('test_execution'))
    nightlySuite.addTest(TALMStaticEvolutionLoadFrictionTestContact('test_execution'))
    nightlySuite.addTest(TALMEvolutionLoadFrictionTestContact('test_execution'))

    # MPC contact test
    nightlySuite.addTest(TBeamAxilSimpleContactTest('test_execution'))
    nightlySuite.addTest(TBeamAxilContactTest('test_execution'))
    nightlySuite.addTest(TBeamAxilTetraContactTest('test_execution'))
    nightlySuite.addTest(TBeamContactTest('test_execution'))
    nightlySuite.addTest(TBeamContactWithTyingTest('test_execution'))
    nightlySuite.addTest(TBeamContactWithFrictionTest('test_execution'))
    nightlySuite.addTest(TPlateTest('test_execution'))

    ### END VALIDATION SUITE ###

    ### BEGIN VALIDATION SUITE ###

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(nightlySuite)

    # ALM frictionless tests
    #nightlySuite.addTest(TALMHertzSphereTestContact('test_execution'))
    validationSuite.addTest(TALMHertzSimpleTestContact('test_execution'))
    validationSuite.addTest(TALMHertzCompleteTestContact('test_execution'))

    # Components ALM frictionless tests
    #nightlySuite.addTest(TComponentsALMHertzSphereTestContact('test_execution'))
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

    # Some large displacement tests
    validationSuite.addTest(TLargeDisplacementPatchTestHexa('test_execution'))
    validationSuite.addTest(TMeshTyingValidationTest('test_execution'))

    # ALM frictionless tests
    validationSuite.addTest(TALMTaylorPatchDynamicTestContact('test_execution'))
    validationSuite.addTest(TALMMeshMovingMatchingTestContact('test_execution'))
    validationSuite.addTest(TALMMeshMovingNotMatchingTestContact('test_execution'))
    #validationSuite.addTest(TALMIroningTestContact('test_execution'))
    #validationSuite.addTest(TALMIroningDieTestContact('test_execution'))
    validationSuite.addTest(TALMLargeDisplacementPatchTestTetra('test_execution'))
    validationSuite.addTest(TALMLargeDisplacementPatchTestHexa('test_execution'))
    validationSuite.addTest(TALMMultiLayerContactTest('test_execution'))
    validationSuite.addTest(TALMSelfContactContactTest('test_execution'))

    # Components ALM frictionless tests
    validationSuite.addTest(TComponentsALMTaylorPatchDynamicTestContact('test_execution'))
    validationSuite.addTest(TComponentsALMMeshMovingMatchingTestContact('test_execution'))
    validationSuite.addTest(TComponentsALMMeshMovingNotMatchingTestContact('test_execution'))
    validationSuite.addTest(TComponentsALMLargeDisplacementPatchTestTetra('test_execution'))
    validationSuite.addTest(TComponentsALMLargeDisplacementPatchTestHexa('test_execution'))
    validationSuite.addTest(TComponentsALMMultiLayerContactTest('test_execution'))
    validationSuite.addTest(TComponentsALMSelfContactContactTest('test_execution'))

    # ALM frictional tests
    validationSuite.addTest(TALMTaylorPatchFrictionalTestContact('test_execution'))
    validationSuite.addTest(TALMMeshMovingMatchingTestFrictionalPureSlipContact('test_execution'))
    validationSuite.addTest(TALMMeshMovingNotMatchingTestFrictionalPureSlipContact('test_execution'))
    validationSuite.addTest(TALMHertzTestFrictionalContact('test_execution'))
    validationSuite.addTest(TALMBlockTestFrictionalContact('test_execution'))

    # MPC contact test
    #validationSuite.addTest(TMultiLayerContactTest('test_execution'))

    ### END VALIDATION ###

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightlySuite) # Already contains the smallSuite
    validationSuite.addTests(allSuite) # Validation contains all

    # Manual list for debugging
    #allSuite.addTests(
        #KratosUnittest.TestLoader().loadTestsFromTestCases([
            #### STANDALONE
            #TTestProcessFactory,
            #TTestDoubleCurvatureIntegration,
            #TTestDynamicSearch,
            #TTestMortarMapperCore,
            #### SMALL
            #TSimplePatchTestTwoDMeshTying,
            #TSimpleSlopePatchTestTwoDMeshTying,
            #TSimplestPatchTestThreeDMeshTying,
            #TSimplestPatchTestThreeDTriQuadMeshTying,
            #TSimplestPatchTestThreeDQuadTriMeshTying,
            #TSimplePatchTestThreeDMeshTying,
            #TALMHyperSimplePatchTestContact,
            #TALMHyperSimplePatchTrianglesTestContact,
            #TALMHyperSimplePatchTestWithEliminationContact,
            #TALMHyperSimplePatchTestWithEliminationWithConstraintContact,
            #TALMHyperSimpleSlopePatchTestContact,
            #TALMTwoDPatchComplexGeomTestContact,
            #TALMTwoDPatchComplexGeomSlopeTestContact,
            #TALMSimplePatchTestContact,
            #TALMSimpleSlopePatchTestContact,
            #TALMSimplePatchNotMatchingATestContact,
            #TALMSimplePatchNotMatchingBTestContact,
            #TALMThreeDSimplestPatchMatchingTestContact,
            #TALMThreeDSimplestPatchTestTriQuadContact,
            #TALMThreeDSimplestPatchTestQuadTriContact,
            #TALMThreeDSimplestPatchMatchingAdaptativeTestContact,
            #TALMThreeDSimplestPatchMatchingSlopeTestContact,
            #TALMThreeDPatchComplexGeomTestContact,
            #TALMTThreeDPatchMatchingTestContact,
            #TALMThreeDPatchNotMatchingTestContact,
            #TPenaltyFrictionlessHyperSimplePatchFrictionalTestContact,
            #TPenaltyThreeDSimplestPatchMatchingTestContact,
            #TComponentsALMHyperSimpleTrianglePatchTestContact,
            #TComponentsALMHyperSimplePatchTestContact,
            #TComponentsALMHyperSimplePatchTestWithEliminationContact,
            #TComponentsALMHyperSimplePatchTestWithEliminationWithConstraintContact,
            #TComponentsALMHyperSimpleSlopePatchTestContact,
            #TComponentsALMTwoDPatchComplexGeomTestContact,
            #TComponentsALMTwoDPatchComplexGeomSlopeTestContact,
            #TComponentsALMSimplePatchTestContact,
            #TComponentsALMSimpleSlopePatchTestContact,
            #TComponentsALMSimplePatchNotMatchingATestContact,
            #TComponentsALMSimplePatchNotMatchingBTestContact,
            #TComponentsALMThreeDSimplestPatchMatchingTestContact,
            #TComponentsALMThreeDSimplestPatchTestTriQuadContact,
            #TComponentsALMThreeDSimplestPatchTestQuadTriContact,
            #TComponentsALMThreeDSimplestPatchMatchingAdaptativeTestContact,
            #TComponentsALMThreeDSimplestPatchMatchingSlopeTestContact,
            #TComponentsALMThreeDPatchComplexGeomTestContact,
            #TComponentsALMTThreeDPatchMatchingTestContact,
            #TComponentsALMThreeDPatchNotMatchingTestContact,
            #TALMHyperSimplePatchFrictionalTestContact,
            #TALMNoFrictionHyperSimplePatchFrictionalTestContact,
            #TALMPerfectStickHyperSimplePatchFrictionalTestContact,
            #TALMThresholdSlipHyperSimplePatchFrictionalTestContact,
            #TALMHyperSimplePatchFrictionalSlipTestContact,
            #TALMHyperSimplePatchFrictionalStickTestContact,
            #TPenaltyNoFrictionHyperSimplePatchFrictionalTestContact,
            #TPenaltyPerfectStickHyperSimplePatchFrictionalTestContact,
            #TPenaltyThresholdSlipHyperSimplePatchFrictionalTestContact,
            #TPenaltyHyperSimplePatchFrictionalSlipTestContact,
            #TPenaltyHyperSimplePatchFrictionalStickTestContact,
            #TTwoDSimplestPatchMatchingTestContact,
            #TTwoDSimplestWithFrictionPatchMatchingTestContact,
            #TThreeDSimplestPatchMatchingTestContact,
            #TThreeDSimplestWithFrictionPatchMatchingTestContact,
            #TThreeDSimplestPatchMatchingSlopeTestContact,
            #TThreeDPatchMatchingTestContact,
            #TThreeDPatchNotMatchingTestContact,
            ##### NIGTHLY
            #TALMTaylorPatchTestContact,
            ######TALMHertzSphereTestContact,  # FIXME: This test requieres the axisymmetric to work (memmory error, correct it)
            #TALMHertzSimpleSphereTestContact,
            #TComponentsALMTaylorPatchTestContact,
            #TALMPureFrictionalTestContact,
            #TALMBasicFrictionTestContact,
            #TALMStaticEvolutionLoadFrictionTestContact,
            #TALMEvolutionLoadFrictionTestContact,
            #TBeamAxilSimpleContactTest,
            #TBeamAxilContactTest,
            #TBeamAxilTetraContactTest,
            #TBeamContactTest,
            #TBeamContactWithTyingTest,
            #TBeamContactWithFrictionTest,
            #TPlateTest,
            ##### VALIDATION
            ######TComponentsALMHertzSphereTestContact,  # FIXME: This test requieres the axisymmetric to work (memmory error, correct it)
            #TALMHertzSimpleTestContact,
            #TALMHertzCompleteTestContact,
            #TComponentsALMHertzSimpleTestContact,
            #TComponentsALMHertzCompleteTestContact,
            #TComponentsALMHertzSimpleSphereTestContact,
            #TExplicitPenaltyThreeDSimplestPatchMatchingTestContact,
            #TALMTaylorPatchDynamicTestContact,
            #TALMMeshMovingMatchingTestContact,
            #TALMMeshMovingNotMatchingTestContact,
            #TALMTaylorPatchFrictionalTestContact,
            #TALMMeshMovingMatchingTestFrictionalPureSlipContact,
            #TALMMeshMovingNotMatchingTestFrictionalPureSlipContact,
            #TALMHertzTestFrictionalContact,
            #TALMBlockTestFrictionalContact,
            #####TALMIroningTestContact,
            #####TALMIroningDieTestContact,
            #TLargeDisplacementPatchTestHexa,
            #TMeshTyingValidationTest,
            #TALMLargeDisplacementPatchTestTetra,
            #TALMLargeDisplacementPatchTestHexa,
            #TALMMultiLayerContactTest,
            #TALMSelfContactContactTest,
            #TComponentsALMTaylorPatchDynamicTestContact,
            #TComponentsALMMeshMovingMatchingTestContact,
            #TComponentsALMMeshMovingNotMatchingTestContact,
            #TComponentsALMLargeDisplacementPatchTestTetra,
            #TComponentsALMLargeDisplacementPatchTestHexa,
            #TComponentsALMMultiLayerContactTest,
            #TComponentsALMSelfContactContactTest,
            ####TMultiLayerContactTest,
        #])
    #)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
