# Definition of the classes for the SMALL TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class undrained_soil_column_2D(TF.TestFactory):
    file_name = "element_tests/undrained_soil_column_2D/undrained_soil_column_2D"
    file_parameters = "element_tests/undrained_soil_column_2D/ProjectParameters.json"


def SetTestSuite(suites):
    small_suite = suites['small']

    small_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            undrained_soil_column_2D
        ])
    )

    return small_suite
