# Importing the Kratos Library
import KratosMultiphysics as KM

if not KM.IsDistributedRun():
    raise Exception("This test script can only be executed in MPI!")

import KratosMultiphysics.KratosUnittest as KratosUnittest

from test_hdf5_model_part_io_mpi import TestCase as TestHDF5ModelPartIO
from test_hdf5_core_mpi import TestOperations as TestHDF5Operations
from test_hdf5_core_mpi import TestFileIO as TestHDF5FileIO
from test_point_set_output_process import TestPointSetOutputProcess
from test_line_output_process import TestLineOutputProcess

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites

    smallSuite = suites['mpi_small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5ModelPartIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5FileIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5Operations]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPointSetOutputProcess, TestLineOutputProcess]))

    nightSuite = suites['mpi_nightly']
    nightSuite.addTests(smallSuite)

    allSuite = suites['mpi_all']
    allSuite.addTests([nightSuite]) # already contains the smallSuite

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
