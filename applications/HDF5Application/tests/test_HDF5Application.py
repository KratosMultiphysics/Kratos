import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_hdf5_model_part_io
from test_hdf5_processes import TestHDF5Processes
import run_cpp_unit_tests

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small']
    smallSuite.addTest(test_hdf5_model_part_io.TestCase('test_HDF5ModelPartIO'))
    smallSuite.addTest(test_hdf5_model_part_io.TestCase('test_HDF5NodalSolutionStepDataIO'))
    smallSuite.addTest(test_hdf5_model_part_io.TestCase('test_HDF5NodalDataValueIO'))
    smallSuite.addTest(test_hdf5_model_part_io.TestCase('test_HDF5ElementDataValueIO'))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5Processes]))

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    allSuite = suites['all']
    allSuite.addTests([nightSuite])

    return suites

if __name__ == '__main__':
    print("Running cpp unit tests ...")
    run_cpp_unit_tests.run()
    print("Finished running cpp unit tests!")
    print("Running python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    print("Finished python tests!")
