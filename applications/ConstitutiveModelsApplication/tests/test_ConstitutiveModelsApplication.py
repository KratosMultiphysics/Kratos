# import Kratos
import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits:

# VALIDATION TESTS
import ValidationTests

from test_modified_cam_clay import TestModifiedCamClayModel

def AssembleTestSuites():

    # Suites to run
    suites = KratosUnittest.KratosSuites

    # SMALL TESTS
    #small_suite = suites['small']

    # NIGHT TESTS
    night_suite = suites['nightly']
    night_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestModifiedCamClayModel]))

    # VALIDATION TESTS
    validation_suite = ValidationTests.SetTestSuite(suites)
    night_suite.addTests(validation_suite)

    # ALL TESTS
    all_suite = suites['all']

    all_suite.addTests(night_suite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
