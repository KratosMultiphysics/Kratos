# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import UsingPyKratos
using_pykratos = UsingPyKratos()

from co_simulation_test_factory import TestSmallCoSimulationCases
from co_simulation_test_factory import TestCoSimulationCases
from test_coupling_interface_data import TestCouplingInterfaceData
from test_data_transfer_operators import TestDataTransferOperators
from test_flower_coupling import TestFLOWerCoupling
from test_coupling_operations import TestScalingOperation
from test_function_callback_utility import TestGenericCallFunction
if not using_pykratos:
    from test_cosim_EMPIRE_API import TestCoSim_EMPIRE_API


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
    ################################################################################
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCouplingInterfaceData]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDataTransferOperators]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestScalingOperation]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGenericCallFunction]))
    if not using_pykratos:
        smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoSim_EMPIRE_API]))


    ################################################################################
    nightSuite = suites['nightly'] # These tests are executed in the nightly build
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSmallCoSimulationCases]))

    nightSuite.addTests(smallSuite)

    ################################################################################
    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoSimulationCases]))
    if not using_pykratos:
        validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestFLOWerCoupling]))

    ################################################################################
    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite
    allSuite.addTests(validationSuite)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
