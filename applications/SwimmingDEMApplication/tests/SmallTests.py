# Definition of the classes for the SMALL TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class candelier(TF.TestFactory):
    file_name = "candelier_tests/CandelierDEM"
    file_parameters = "candelier_tests/ProjectParametersDEM.json"

def SetTestSuite(suites):
    small_suite = suites['small']

    small_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            candelier
        ])
    )

    return small_suite
