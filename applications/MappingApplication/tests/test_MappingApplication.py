# import Kratos
import KratosMultiphysics
import KratosMultiphysics.MappingApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suits
import test_nearest_neighbor_mapper
import test_nearest_neighbor_mapper_iga
import test_nearest_element_mapper
import test_barycentric_mapper
import test_projection_3d_2d_mapper
import test_coupling_geometry_mapper

from test_patch_test_mappers import TestPatchTestMappers

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

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
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_coupling_geometry_mapper.TestIgaFEMCouplingGeometryMapper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_coupling_geometry_mapper.TestDualMortarIgaFEMCouplingGeometryMapper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_coupling_geometry_mapper.TestDualMortarCouplingGeometryMapper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_coupling_geometry_mapper.TestSlaveOriginCouplingGeometryMapper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_coupling_geometry_mapper.TestComputeMappingMatrixCouplingGeometryMapper]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsLine]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsVolume]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BladeMapping]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.CornerCaseNearestNeighbor]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper_iga.BasicTestsLineMappingIGAFEM]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper_iga.BasicTestsSurfaceMappingIGAFEM]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsLine]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsVolume]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BladeMapping]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsLineMappingIGAFEM]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsSurfaceMappingIGAFEM]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BasicTestsLine]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BasicTestsVolume]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BladeMapping]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_projection_3d_2d_mapper.Projection3D2DMapperNearestNeighborOrigin2D]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_projection_3d_2d_mapper.Projection3D2DMapperNearestElementOrigin2D]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_projection_3d_2d_mapper.Projection3D2DMapperNearestNeighborDestination2D]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_projection_3d_2d_mapper.Projection3D2DMapperNearestElementDestination2D]))

    # Create a test suit with the selected tests
    # nightSuite will contain the following tests:
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsLineSwitchedSides]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsSurface]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsSurfaceSwitchedSides]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsVolumeSwitchedSides]))
    
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper_iga.BasicTestsLineMappingIGAFEM]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper_iga.BasicTestsSurfaceMappingIGAFEM]))

    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsLineSwitchedSides]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsSurface]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsSurfaceSwitchedSides]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsVolumeSwitchedSides]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsLineMappingIGAFEM]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsSurfaceMappingIGAFEM]))

    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BasicTestsLineSwitchedSides]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BasicTestsSurface]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BasicTestsSurfaceSwitchedSides]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BasicTestsVolumeSwitchedSides]))

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
