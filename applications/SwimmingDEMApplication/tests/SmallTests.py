# Definition of the classes for the SMALL TESTS

#Iimport Kratos
import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.SwimmingDEMApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def SetTestSuite(suites):
    small_suite = suites['small']

    small_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
        ])
    )

    return small_suite

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites
    night_suite = SetTestSuite(suites)
    suites['all'].addTests(night_suite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.DETAIL)
    KratosUnittest.runTests(AssembleTestSuites())

