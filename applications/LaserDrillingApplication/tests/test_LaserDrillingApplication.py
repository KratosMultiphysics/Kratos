# import Kratos
import KratosMultiphysics as KM

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

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

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']

    # Create a test suit with the selected tests plus all small tests
    nightlySuite = suites['nightly']
    nightlySuite.addTests(smallSuite)

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightlySuite)

    return suites

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
