# Importing the Kratos Library
import KratosMultiphysics as KM

if not KM.IsDistributedRun():
    raise Exception("This test script can only be executed in MPI!")

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suits
import test_model_part_utils
import test_container_expression_utils
import test_container_expression
import test_collective_expressions
import process_tests.test_optimization_problem_ascii_output_process

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

    ### Small MPI tests ########################################################
    smallMPISuite = suites['mpi_small']

    # adding custom process tests
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model_part_utils.TestOptAppModelPartUtils]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression_utils.TestContainerExpressionUtils]))

    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression.TestConditionPropertiesExpression]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression.TestElementPropertiesExpression]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_collective_expressions.TestCollectiveExpressions]))

    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([process_tests.test_optimization_problem_ascii_output_process.TestOptimizationProblemAsciiOutputProcess]))

    ### Nightly MPI tests ######################################################
    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTests(smallMPISuite)

    ### Full MPI set ###########################################################
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite) # already contains the smallMPISuite

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
