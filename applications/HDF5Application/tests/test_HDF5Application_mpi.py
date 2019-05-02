import KratosMultiphysics.KratosUnittest as KratosUnittest

from test_hdf5_model_part_io_mpi import TestCase as TestHDF5ModelPartIO

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites

    smallSuite = suites['mpi_small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5ModelPartIO]))

    nightSuite = suites['mpi_nightly']
    nightSuite.addTests(smallSuite)

    allSuite = suites['mpi_all']
    allSuite.addTests([nightSuite]) # already contains the smallSuite

    # temporary until mpi-testing is properly implemented
    auxallSuite = suites['all']
    auxallSuite.addTests(allSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
