# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.FSIApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from KratosExecuteConvergenceAcceleratorSpringMPITest import KratosExecuteConvergenceAcceleratorSpringMPITest as TConvergenceAcceleratorSpringTest

## NIGTHLY TESTS

## VALIDATION TESTS
from ValidationTests import MokBenchmarkTest as TMokBenchmark

def AssambleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    ### Small MPI tests ########################################################
    smallMPISuite = suites['mpi_small']
    smallMPISuite.addTest(TConvergenceAcceleratorSpringTest('test_aitken_accelerator_constant_forces'))
    smallMPISuite.addTest(TConvergenceAcceleratorSpringTest('test_aitken_accelerator_variable_stiffness'))
    smallMPISuite.addTest(TConvergenceAcceleratorSpringTest('test_aitken_accelerator_ghost_nodes'))

    ### Nightly MPI tests ######################################################
    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTests(smallMPISuite)

    ### Full MPI set ###########################################################
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite)

    allSuite = suites['all']
    allSuite.addTests(allMPISuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests( AssambleTestSuites() )
