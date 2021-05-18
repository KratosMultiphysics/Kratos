# Definition of the classes for the NIGHTLY TESTS

#Iimport Kratos
import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.SwimmingDEMApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import SmallTests

# List of tests that are available
available_tests = []

def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(available_tests))

    return night_suite

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites
    small_suite = SmallTests.SetTestSuite(suites)
    suites['all'].addTests(small_suite)
    night_suite = SetTestSuite(suites)
    suites['all'].addTests(night_suite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
