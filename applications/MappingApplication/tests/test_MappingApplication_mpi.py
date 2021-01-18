# import Kratos
import KratosMultiphysics
if not KratosMultiphysics.IsDistributedRun():
    raise Exception("These tests can only be executed in MPI / distributed!")

import KratosMultiphysics.MappingApplication

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

    ### Nightly MPI tests ######################################################
    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTests(smallMPISuite)
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsLine]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsLineSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsSurface]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsSurfaceSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsVolume]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBasicTestsVolumeSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.NearestNeighborBladeMapping]))

    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsLine]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsLineSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsSurface]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsSurfaceSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsVolume]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBasicTestsVolumeSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.NearestElementBladeMapping]))

    ### Full MPI set ###########################################################
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite) # already contains the smallMPISuite

    # can be removed after the cmd-line of the testing accepts "--using-mpi"
    allSuite = suites['all']
    allSuite.addTests(allMPISuite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    # TODO throw if --using-mpi is not being passed!
    KratosUnittest.runTests(AssembleTestSuites())
