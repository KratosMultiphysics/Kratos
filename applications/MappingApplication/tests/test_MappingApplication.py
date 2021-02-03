# import Kratos
import KratosMultiphysics
import KratosMultiphysics.MappingApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

import run_cpp_unit_tests

# Import the tests or test_classes to create the suits
import test_nearest_neighbor_mapper
import test_nearest_element_mapper
import test_coupling_geometry_mapper

from test_patch_test_mappers import TestPatchTestMappers

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

    # Create a test suit with the selected tests (Small tests):
    # smallSuite will contain the following tests:
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPatchTestMappers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_coupling_geometry_mapper.TestCouplingGeometryMapper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_coupling_geometry_mapper.TestDualMortarCouplingGeometryMapper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_coupling_geometry_mapper.TestSlaveOriginCouplingGeometryMapper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_coupling_geometry_mapper.TestComputeMappingMatrixCouplingGeometryMapper]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsLine]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsSurface]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsVolume]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBladeMapping]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsLine]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsSurface]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsVolume]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBladeMapping]))

    # Create a test suit with the selected tests
    # nightSuite will contain the following tests:
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsLineSwitchedSides]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsSurfaceSwitchedSides]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsVolumeSwitchedSides]))

    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsLineSwitchedSides]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsSurfaceSwitchedSides]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsVolumeSwitchedSides]))

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    print("Running cpp unit tests for MappingApplication...")
    run_cpp_unit_tests.run()
    print("Finished running cpp unit tests!")
    print("Running python tests for MappingApplication...")
    KratosUnittest.runTests(AssembleTestSuites())
    print("Finished python tests!")