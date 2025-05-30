# import Kratos
import KratosMultiphysics
if not KratosMultiphysics.IsDistributedRun():
    raise Exception("This test script can only be executed in MPI!")

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suits
import test_nearest_neighbor_mapper
import test_nearest_element_mapper
import test_barycentric_mapper
import test_projection_3d_2d_mapper

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

    ### Small MPI tests ########################################################
    smallMPISuite = suites['mpi_small']
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsLine]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsVolume]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BladeMappingSerialModelPart]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BladeMappingUnevenRanks]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.CornerCaseNearestNeighbor]))

    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsLine]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsVolume]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BladeMappingSerialModelPart]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BladeMappingUnevenRanks]))

    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BasicTestsLine]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BasicTestsVolume]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BladeMappingSerialModelPart]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BladeMappingUnevenRanks]))

    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_projection_3d_2d_mapper.Projection3D2DMapperNearestNeighborOrigin2D]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_projection_3d_2d_mapper.Projection3D2DMapperNearestElementOrigin2D]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_projection_3d_2d_mapper.Projection3D2DMapperNearestNeighborDestination2D]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_projection_3d_2d_mapper.Projection3D2DMapperNearestElementDestination2D]))

    ### Nightly MPI tests ######################################################
    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTests(smallMPISuite)
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsLineSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsSurface]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsSurfaceSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsVolumeSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BladeMapping]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BladeMappingAllRanksExceptLast]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BladeMappingAllRanksExceptFirst]))

    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsLineSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsSurface]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsSurfaceSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsVolumeSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BladeMapping]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BladeMappingAllRanksExceptLast]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BladeMappingAllRanksExceptFirst]))

    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BasicTestsLineSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BasicTestsSurface]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BasicTestsSurfaceSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BasicTestsVolumeSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BladeMapping]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BladeMappingAllRanksExceptLast]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_barycentric_mapper.BladeMappingAllRanksExceptFirst]))

    ### Full MPI set ###########################################################
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite) # already contains the smallMPISuite

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
