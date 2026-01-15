# import Kratos
import KratosMultiphysics

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suites
import test_vtu_output
import test_sparse_matrix_linear_operator

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

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sparse_matrix_linear_operator.TestSparseMatrixLinearOperator]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_vtu_output.TestVtuOutput2D]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_vtu_output.TestVtuOutput3D]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_vtu_output.TestVtuOutput]))

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
