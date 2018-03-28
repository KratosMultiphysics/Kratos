# Definition of the classes for the VALIDATION TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

# to set

def SetTestSuite(suites):
    validation_suite = suites['validation']

    validation_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
        ])
    )

    return validation_suite
