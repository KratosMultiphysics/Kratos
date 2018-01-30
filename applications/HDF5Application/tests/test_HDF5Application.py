import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_hdf5_model_part_io

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small']
    smallSuite.addTest(test_hdf5_model_part_io.TestCase('test_HDF5ModelPartIO'))
    smallSuite.addTest(test_hdf5_model_part_io.TestCase('test_HDF5NodalSolutionStepDataIO'))

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    allSuite = suites['all']
    allSuite.addTests([smallSuite])

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
