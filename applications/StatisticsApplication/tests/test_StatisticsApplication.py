# import Kratos
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.StatisticsApplication

import os, subprocess

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
import run_cpp_unit_tests
from test_norms import NormTests
from test_spatial_methods import SpatialMethodTests
from test_temporal_sum_method import TemporalSumMethodTests
from test_temporal_mean_method import TemporalMeanMethodTests
from test_temporal_variance_method import TemporalVarianceMethodTests
from test_temporal_min_method import TemporalMinMethodTests
from test_temporal_max_method import TemporalMaxMethodTests
from test_temporal_rms_method import TemporalRootMeanSquareMethodTests
from test_spatial_statistics_process import SpatialStatisticsProcessTest

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
    # - testSmallExample
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([NormTests]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([SpatialMethodTests]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TemporalSumMethodTests]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TemporalMeanMethodTests]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TemporalVarianceMethodTests]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TemporalMinMethodTests]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TemporalMaxMethodTests]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TemporalRootMeanSquareMethodTests]))

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([SpatialStatisticsProcessTest]))

    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    if kratos_utilities.IsMPIAvailable() and kratos_utilities.CheckIfApplicationsAvailable("MetisApplication"):
        KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning mpi python tests ...")
        p = subprocess.Popen(
            ["mpiexec", "-np", "2", "python3", "test_StatisticsApplication_mpi.py", "--using-mpi"],
            stdout=subprocess.PIPE,
            cwd=os.path.dirname(os.path.abspath(__file__)))
        p.wait()
        KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished mpi python tests!")
    else:
        KratosMultiphysics.Logger.PrintInfo("Unittests", "\nSkipping mpi python tests due to missing dependencies")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
