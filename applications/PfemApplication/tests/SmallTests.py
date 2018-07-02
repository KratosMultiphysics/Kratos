# Definition of the classes for the SMALL TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class newtonian_sloshing_2D(TF.TestFactory):
    file_name = "fluid_tests/newtonian/sloshing_2D"
    file_parameters = None


def SetTestSuite(suites):
    small_suite = suites['small']

    small_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            newtonian_sloshing_2D
        ])
    )

    return small_suite
