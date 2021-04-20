# import Kratos
import KratosMultiphysics
import KratosMultiphysics.ConstitutiveLawsApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from test_perfect_plasticity_implementation_verification import TestPerfectPlasticityImplementationVerification

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

    # Create a test suit with the selected tests (Nightly tests):
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPerfectPlasticityImplementationVerification]))


    ### Adding Validation Tests
    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
