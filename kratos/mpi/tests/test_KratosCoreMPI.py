from __future__ import print_function, absolute_import, division

# import Kratos
import KratosMultiphysics
import KratosMultiphysics.mpi as MPI #TODO: do not import the so directly (but I need a nice Python module first)

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suites
import test_mpi_communicator_set_up
import test_mpi_data_communicator_python

def AssembleMPITestSuites():
    ''' Populates the test suites to run.

    Populates the MPI test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mpi_communicator_set_up.TestMPICommunicatorSetUp]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mpi_data_communicator_python.TestMPIDataCommunicatorPython]))

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleMPITestSuites())
