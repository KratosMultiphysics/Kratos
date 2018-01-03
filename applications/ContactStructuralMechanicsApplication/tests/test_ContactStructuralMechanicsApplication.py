# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ContactStructuralMechanicsApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
# Exact integration tests
from test_double_curvature_integration import TestDoubleCurvatureIntegration as TTestDoubleCurvatureIntegration
from test_mortar_mapper import TestMortarMapping as TTestMortarMapping

# Mesh tying tests
from SmallTests import SimplePatchTestTwoDMeshTying      as TSimplePatchTestTwoDMeshTying
from SmallTests import SimpleSlopePatchTestTwoDMeshTying as TSimpleSlopePatchTestTwoDMeshTying
from SmallTests import SimplestPatchTestThreeDMeshTying  as TSimplestPatchTestThreeDMeshTying
from SmallTests import SimplePatchTestThreeDMeshTying    as TSimplePatchTestThreeDMeshTying

# ALM frictionless tests
from SmallTests import ALMHyperSimplePatchTestContact                      as TALMHyperSimplePatchTestContact
from SmallTests import ALMHyperSimpleSlopePatchTestContact                 as TALMHyperSimpleSlopePatchTestContact
from SmallTests import ALMTwoDPatchComplexGeomTestContact                  as TALMTwoDPatchComplexGeomTestContact
from SmallTests import ALMTwoDPatchComplexGeomSlopeTestContact             as TALMTwoDPatchComplexGeomSlopeTestContact
from SmallTests import ALMSimplePatchTestContact                           as TALMSimplePatchTestContact
from SmallTests import ALMSimpleSlopePatchTestContact                      as TALMSimpleSlopePatchTestContact
from SmallTests import ALMSimplePatchNotMatchingATestContact               as TALMSimplePatchNotMatchingATestContact
from SmallTests import ALMSimplePatchNotMatchingBTestContact               as TALMSimplePatchNotMatchingBTestContact
from SmallTests import ALMThreeDSimplestPatchMatchingTestContact           as TALMThreeDSimplestPatchMatchingTestContact
from SmallTests import ALMThreeDSimplestPatchMatchingAdaptativeTestContact as TALMThreeDSimplestPatchMatchingAdaptativeTestContact
from SmallTests import ALMThreeDSimplestPatchMatchingSlopeTestContact      as TALMThreeDSimplestPatchMatchingSlopeTestContact
from SmallTests import ALMThreeDPatchComplexGeomTestContact                as TALMThreeDPatchComplexGeomTestContact
from SmallTests import ALMThreeDPatchMatchingTestContact                   as TALMTThreeDPatchMatchingTestContact
from SmallTests import ALMThreeDPatchNotMatchingTestContact                as TALMThreeDPatchNotMatchingTestContact

## NIGTHLY TESTS
# ALM frictionless tests
from NightlyTests import ALMMeshMovingMatchingTestContact    as TALMMeshMovingMatchingTestContact
from NightlyTests import ALMMeshMovingNotMatchingTestContact as TALMMeshMovingNotMatchingTestContact
from NightlyTests import ALMTaylorPatchTestContact           as TALMTaylorPatchTestContact
from NightlyTests import ALMTaylorPatchDynamicTestContact    as TALMTaylorPatchDynamicTestContact
from NightlyTests import ALMHertzSimpleTestContact           as TALMHertzSimpleTestContact
from NightlyTests import ALMHertzSimpleSphereTestContact     as TALMHertzSimpleSphereTestContact
from NightlyTests import ALMHertzSphereTestContact           as TALMHertzSphereTestContact
from NightlyTests import ALMHertzCompleteTestContact         as TALMHertzCompleteTestContact

## VALIDATION TESTS
# ALM frictionless tests
from ValidationTests import ALMIroningTestContact    as TALMIroningTestContact
from ValidationTests import ALMIroningDieTestContact as TALMIroningDieTestContact
from ValidationTests import LargeDisplacementPatchTestHexa as TLargeDisplacementPatchTestHexa
from ValidationTests import ALMLargeDisplacementPatchTestTetra as TALMLargeDisplacementPatchTestTetra
from ValidationTests import ALMLargeDisplacementPatchTestHexa as TALMLargeDisplacementPatchTestHexa

def AssambleTestSuites():
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
    # Exact integration tests
    smallSuite.addTest(TTestDoubleCurvatureIntegration('test_double_curvature_integration_triangle'))
    smallSuite.addTest(TTestDoubleCurvatureIntegration('test_double_curvature_integration_quad'))
    smallSuite.addTest(TTestDoubleCurvatureIntegration('test_moving_mesh_integration_quad'))
    
    # Mortar mapping
    smallSuite.addTest(TTestMortarMapping('test_basic_mortar_mapping_triangle'))
    smallSuite.addTest(TTestMortarMapping('test_basic_mortar_mapping_quad'))
    smallSuite.addTest(TTestMortarMapping('test_less_basic_mortar_mapping_triangle'))
    smallSuite.addTest(TTestMortarMapping('test_less_basic_2_mortar_mapping_triangle'))
    smallSuite.addTest(TTestMortarMapping('test_simple_curvature_mortar_mapping_triangle'))
    smallSuite.addTest(TTestMortarMapping('test_mortar_mapping_triangle'))
    smallSuite.addTest(TTestMortarMapping('test_mortar_mapping_quad'))
    
    # Mesh tying tests 
    smallSuite.addTest(TSimplePatchTestTwoDMeshTying('test_execution'))
    smallSuite.addTest(TSimpleSlopePatchTestTwoDMeshTying('test_execution'))
    smallSuite.addTest(TSimplestPatchTestThreeDMeshTying('test_execution'))
    smallSuite.addTest(TSimplePatchTestThreeDMeshTying('test_execution'))
    
    # ALM frictionless tests
    smallSuite.addTest(TALMHyperSimplePatchTestContact('test_execution'))
    smallSuite.addTest(TALMHyperSimpleSlopePatchTestContact('test_execution'))
    smallSuite.addTest(TALMTwoDPatchComplexGeomTestContact('test_execution'))
    smallSuite.addTest(TALMTwoDPatchComplexGeomSlopeTestContact('test_execution'))
    smallSuite.addTest(TALMSimplePatchTestContact('test_execution'))
    smallSuite.addTest(TALMSimpleSlopePatchTestContact('test_execution'))
    smallSuite.addTest(TALMSimplePatchNotMatchingATestContact('test_execution'))
    smallSuite.addTest(TALMSimplePatchNotMatchingBTestContact('test_execution'))
    smallSuite.addTest(TALMThreeDSimplestPatchMatchingTestContact('test_execution'))
    smallSuite.addTest(TALMThreeDSimplestPatchMatchingAdaptativeTestContact('test_execution'))
    smallSuite.addTest(TALMThreeDSimplestPatchMatchingSlopeTestContact('test_execution'))
    smallSuite.addTest(TALMThreeDPatchComplexGeomTestContact('test_execution'))
    smallSuite.addTest(TALMTThreeDPatchMatchingTestContact('test_execution'))
    smallSuite.addTest(TALMThreeDPatchNotMatchingTestContact('test_execution'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTest(TALMMeshMovingMatchingTestContact('test_execution'))
    nightSuite.addTest(TALMMeshMovingNotMatchingTestContact('test_execution'))
    nightSuite.addTest(TALMTaylorPatchTestContact('test_execution'))
    nightSuite.addTest(TALMTaylorPatchDynamicTestContact('test_execution'))
    nightSuite.addTest(TALMHertzSimpleSphereTestContact('test_execution'))
    nightSuite.addTest(TALMHertzSphereTestContact('test_execution'))
    nightSuite.addTest(TALMHertzSimpleTestContact('test_execution'))
    nightSuite.addTest(TALMHertzCompleteTestContact('test_execution'))
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']
    validationSuite.addTests(nightSuite)
    #validationSuite.addTest(TALMIroningTestContact('test_execution'))
    #validationSuite.addTest(TALMIroningDieTestContact('test_execution'))
    validationSuite.addTest(TLargeDisplacementPatchTestHexa('test_execution'))
    validationSuite.addTest(TALMLargeDisplacementPatchTestTetra('test_execution'))
    validationSuite.addTest(TALMLargeDisplacementPatchTestHexa('test_execution'))
    
    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            ## SMALL
            TTestDoubleCurvatureIntegration,
            TTestMortarMapping,
            TSimplePatchTestTwoDMeshTying,
            TSimpleSlopePatchTestTwoDMeshTying,
            TSimplestPatchTestThreeDMeshTying,
            TSimplePatchTestThreeDMeshTying,
            TALMHyperSimplePatchTestContact,
            TALMHyperSimpleSlopePatchTestContact,
            TALMTwoDPatchComplexGeomTestContact,
            TALMTwoDPatchComplexGeomSlopeTestContact,
            TALMSimplePatchTestContact,
            TALMSimpleSlopePatchTestContact, 
            TALMSimplePatchNotMatchingATestContact,
            TALMSimplePatchNotMatchingBTestContact,
            TALMThreeDSimplestPatchMatchingTestContact,
            TALMThreeDSimplestPatchMatchingAdaptativeTestContact,
            TALMThreeDSimplestPatchMatchingSlopeTestContact,
            TALMThreeDPatchComplexGeomTestContact,
            TALMTThreeDPatchMatchingTestContact,
            TALMThreeDPatchNotMatchingTestContact,
            ## NIGTHLY
            TALMMeshMovingMatchingTestContact,
            TALMMeshMovingNotMatchingTestContact,
            TALMTaylorPatchTestContact,
            TALMTaylorPatchDynamicTestContact,
            TALMHertzSimpleTestContact,
            TALMHertzSimpleSphereTestContact,
            ##TALMHertzSphereTestContact,  # FIXME: This test requieres the axisymmetric to work (memmory error, correct it)
            TALMHertzCompleteTestContact,
            ## VALIDATION
            ##TALMIroningTestContact,
            ##TALMIroningDieTestContact,
            TLargeDisplacementPatchTestHexa,
            TALMLargeDisplacementPatchTestTetra,
            TALMLargeDisplacementPatchTestHexa,
        ])
    )

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
