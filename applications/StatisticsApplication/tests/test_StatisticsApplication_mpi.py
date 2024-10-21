# import Kratos
import KratosMultiphysics
import KratosMultiphysics.StatisticsApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suits

# Shell tests
from test_spatial_methods import SpatialMethodTests
from test_spatial_statistics_process import SpatialStatisticsProcessTest
from test_temporal_sum_method import TemporalSumMethodTests
from test_temporal_mean_method import TemporalMeanMethodTests
from test_temporal_variance_method import TemporalVarianceMethodTests
from test_temporal_min_method import TemporalMinMethodTests
from test_temporal_max_method import TemporalMaxMethodTests
from test_temporal_rms_method import TemporalRootMeanSquareMethodTests

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
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([SpatialMethodTests]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([SpatialStatisticsProcessTest]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TemporalSumMethodTests]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TemporalMeanMethodTests]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TemporalVarianceMethodTests]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TemporalMinMethodTests]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TemporalMaxMethodTests]))
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TemporalRootMeanSquareMethodTests]))

    ### Nightly MPI tests ######################################################
    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTests(smallMPISuite)

    ### Full MPI set ###########################################################
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite) # already contains the smallMPISuite

    # can be removed after the cmd-line of the testing accepts "--using-mpi"
    allSuite = suites['all']
    allSuite.addTests(allMPISuite)

    return suites


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
