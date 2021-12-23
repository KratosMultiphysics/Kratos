# import Kratos
import KratosMultiphysics
if not KratosMultiphysics.IsDistributedRun():
    raise Exception("This test script can only be executed in MPI!")

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suits
import run_cpp_mpi_unit_tests
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
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsLine]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsSurface]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsVolume]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BladeMappingSerialModelPart]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BladeMappingUnevenRanks]))

    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsLine]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsSurface]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsVolume]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BladeMappingSerialModelPart]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BladeMappingUnevenRanks]))

    ### Nightly MPI tests ######################################################
    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTests(smallMPISuite)
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsLineSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsSurfaceSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BasicTestsVolumeSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BladeMapping]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BladeMappingAllRanksExceptLast]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_neighbor_mapper.BladeMappingAllRanksExceptFirst]))

    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsLineSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsSurfaceSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BasicTestsVolumeSwitchedSides]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BladeMapping]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BladeMappingAllRanksExceptLast]))
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_nearest_element_mapper.BladeMappingAllRanksExceptFirst]))

    ### Full MPI set ###########################################################
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite) # already contains the smallMPISuite

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    run_cpp_mpi_unit_tests.run() # Application-MPI tests are currently not run automatically, hence this hack is temporarily required
    KratosUnittest.runTests(AssembleTestSuites())
