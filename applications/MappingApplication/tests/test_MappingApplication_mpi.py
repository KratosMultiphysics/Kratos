# import Kratos
import KratosMultiphysics
if not KratosMultiphysics.IsDistributedRun():
    raise Exception("This test script can only be executed in MPI!")

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suits
import test_nearest_neighbor_mapper
import test_nearest_element_mapper

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

    ### Small MPI tests ########################################################
    smallMPISuite = suites['mpi_small']
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsLine]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsSurface]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsVolume]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBladeMapping]))

    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsLine]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsSurface]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsVolume]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBladeMapping]))

    ### Nightly MPI tests ######################################################
    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTests(smallMPISuite)
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsLineSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsSurfaceSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsVolumeSwitchedSides]))

    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsLineSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsSurfaceSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsVolumeSwitchedSides]))

    ### Full MPI set ###########################################################
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite) # already contains the smallMPISuite

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
