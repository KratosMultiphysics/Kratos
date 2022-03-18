import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.ThermalDEMApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import SmallTests
import NightTests

KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

def AssembleTestSuites():
    # Suites to run
    suites = KratosUnittest.KratosSuites
    
    # Small and nightly tests
    small_suite = SmallTests.SetTestSuite(suites)
    night_suite = NightTests.SetTestSuite(suites)
    
    # Include small suite in night suite
    night_suite.addTests(small_suite)
    
    # All tests
    all_suite = suites['all']
    all_suite.addTests(night_suite)
    
    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
