# Importing the Kratos Library
import KratosMultiphysics as KM

if not KM.IsDistributedRun():
    raise Exception("This test script can only be executed in MPI!")

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suits
from test_model_part_utils import TestModelPartUtils
from test_container_expression_utils import TestContainerExpressionUtils
from test_container_expression import TestConditionPropertiesExpression
from test_container_expression import TestElementPropertiesExpression

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
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestModelPartUtils]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestContainerExpressionUtils]))

    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConditionPropertiesExpression]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestElementPropertiesExpression]))

    ### Nightly MPI tests ######################################################
    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTests(smallMPISuite)

    ### Full MPI set ###########################################################
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite) # already contains the smallMPISuite

    return suites


if __name__ == '__main__':
    KM.Tester.SetVerbosity(KM.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    KratosUnittest.runTests(AssembleTestSuites())
