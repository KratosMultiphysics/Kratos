# import Kratos
import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.SwimmingDEMApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits:
import SmallTests
import NightTests
KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

def AssembleTestSuites():

    # Suites to run
    suites = KratosUnittest.KratosSuites

    # SMALL TESTS
    small_suite = SmallTests.SetTestSuite(suites)

    # NIGHTLY TESTS
    night_suite = NightTests.SetTestSuite(suites)

    # include small suite in night suite
    night_suite.addTests(small_suite)

    # ALL TESTS
    all_suite = suites['all']

    all_suite.addTests(night_suite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.DETAIL)
    KratosUnittest.runTests(AssembleTestSuites())


