# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.FSIApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from convergence_accelerator_spring_MPI_test import ConvergenceAcceleratorSpringMPITest

## NIGTHLY TESTS

## VALIDATION TESTS

def AssembleTestSuites():
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
    smallMPISuite.addTest(ConvergenceAcceleratorSpringMPITest('test_aitken_accelerator_constant_forces'))
    smallMPISuite.addTest(ConvergenceAcceleratorSpringMPITest('test_aitken_accelerator_variable_stiffness'))
    smallMPISuite.addTest(ConvergenceAcceleratorSpringMPITest('test_aitken_accelerator_ghost_nodes'))
    smallMPISuite.addTest(ConvergenceAcceleratorSpringMPITest('test_mvqn_recursive_accelerator_constant_forces'))
    smallMPISuite.addTest(ConvergenceAcceleratorSpringMPITest('test_mvqn_recursive_accelerator_variable_stiffness'))
    smallMPISuite.addTest(ConvergenceAcceleratorSpringMPITest('test_mvqn_recursive_accelerator_ghost_nodes'))

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
    KratosUnittest.runTests( AssembleTestSuites() )
