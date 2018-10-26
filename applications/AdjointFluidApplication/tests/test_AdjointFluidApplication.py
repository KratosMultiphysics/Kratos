import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

## NIGTHLY TESTS

## VALIDATION TESTS

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['small']


    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nightly and you can use to validate
    validationSuite = suites['validation']

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests([smallSuite])

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
