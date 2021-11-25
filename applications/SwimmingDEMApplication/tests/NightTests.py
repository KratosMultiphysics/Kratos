# Definition of the classes for the NIGHTLY TESTS

#Iimport Kratos
import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.SwimmingDEMApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import SmallTests
import GentleInjectionTestFactory as GentleTF

class gentle_injection_test(GentleTF.GentleInjectionTestFactory):
     file_name = "gentle_injection_tests/cube_cavity_with_inlet"
     file_parameters = "gentle_injection_tests/ProjectParameters.json"

# List of tests that are available
available_tests = []
available_tests += [test_class for test_class in GentleTF.GentleInjectionTestFactory.__subclasses__()]

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
