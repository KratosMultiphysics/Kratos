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
    small_suite = suites['small']
    night_suite = suites['nightly']
    all_suite = suites['all']
    validation_suite = suites['validation']

    small_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestModifiedCamClayModel]))
    night_suite = suites['small']
    all_suite = suites['small']
    validation_suite = suites['small']

    # NIGHT TESTS

    # ALL TESTS


    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")


