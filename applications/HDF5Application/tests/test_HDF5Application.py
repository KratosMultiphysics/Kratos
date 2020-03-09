import subprocess
from os.path import dirname
from os.path import abspath


import KratosMultiphysics.KratosUnittest as KratosUnittest
from test_hdf5_core import TestFileIO as TestHDF5FileIO
from test_hdf5_core import TestOperations as TestHDF5Operations
from test_hdf5_core import TestControllers as TestHDF5Controllers
from test_hdf5_model_part_io import TestCase as TestHDF5ModelPartIO
from test_hdf5_processes import TestHDF5Processes
from test_hdf5_xdmf import TestTryOpenH5File
from test_hdf5_xdmf import TestCreateXdmfSpatialGrid
from test_hdf5_xdmf import TestXdmfNodalResults
from test_hdf5_xdmf import TestXdmfElementResults
from test_hdf5_xdmf import TestXdmfConditionResults
from test_hdf5_xdmf import TestXdmfResults
from test_hdf5_xdmf import TestTimeLabel
from test_hdf5_xdmf import TestFindMatchingFiles
from test_hdf5_xdmf import TestCreateXdmfTemporalGridFromMultifile


def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5FileIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5Operations]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5Controllers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5ModelPartIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5Processes]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTryOpenH5File]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCreateXdmfSpatialGrid]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestXdmfNodalResults]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestXdmfElementResults]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestXdmfConditionResults]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestXdmfResults]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTimeLabel]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestFindMatchingFiles]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCreateXdmfTemporalGridFromMultifile]))
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    allSuite = suites['all']
    allSuite.addTests([nightSuite])
    return suites


def run_cpp_unit_tests():
    """Run the c++ unit tests.

    This runs a subprocess because the HDF5 used by h5py clashed with the HDF5
    linked to the c++ unit tests on some systems. h5py is used by some XDMF tests.
    """
    # We set cwd in case the script is run from another directory. This is needed
    # when testing from the core.
    out_bytes = subprocess.check_output(
        ['python3', 'run_cpp_unit_tests.py'], cwd=abspath(dirname(__file__)))
    return out_bytes.decode('utf-8')


if __name__ == '__main__':
    print("\nRunning cpp unit tests ...")
    cpp_test_text = run_cpp_unit_tests()
    print(cpp_test_text)
    print("Finished running cpp unit tests!")
    print("\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    print("Finished python tests!")
