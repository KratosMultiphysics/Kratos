import subprocess
from os.path import dirname
from os.path import abspath


import KratosMultiphysics.KratosUnittest as KratosUnittest
from test_hdf5_core import TestFileIO as TestHDF5FileIO
from test_hdf5_core import TestOperations as TestHDF5Operations
from test_hdf5_core import TestControllers as TestHDF5Controllers
from test_hdf5_model_part_io import TestCase as TestHDF5ModelPartIO
from test_hdf5_processes import TestHDF5Processes
from test_hdf5_xdmf import TestCreateXdmfSpatialGrid
from test_hdf5_xdmf import TestXdmfNodalResults
from test_hdf5_xdmf import TestXdmfElementResults
from test_hdf5_xdmf import TestXdmfConditionResults
from test_vertex import TestVertex
from test_point_set_output_process import TestPointSetOutputProcess as TestHDF5PointSetOutputProcess
from test_line_output_process import TestLineOutputProcess as TestHDF5LineOutputProcess
from test_pattern import TestGetMachingEntitiesString
from test_expression_io import TestExpressionIO
from test_hdf5_mesh_location_container import TestMeshLocationContainer
from test_dataset_generator import TestDatasetGenerator

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5FileIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5Operations]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5Controllers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5ModelPartIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHDF5Processes]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCreateXdmfSpatialGrid]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestXdmfNodalResults]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestXdmfElementResults]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestXdmfConditionResults]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestVertex, TestHDF5PointSetOutputProcess, TestHDF5LineOutputProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestExpressionIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGetMachingEntitiesString]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMeshLocationContainer]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDatasetGenerator]))
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
