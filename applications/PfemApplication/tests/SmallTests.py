# Definition of the classes for the SMALL TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class newtonian_dam_break_2D(TF.TestFactory):
    file_name = "fluid_tests/newtonian/dam_break_2D"
    file_parameters = None

def SetTestSuite(suites):
    small_suite = suites['small']

    small_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            newtonian_dam_break_2D
        ])
    )

    return small_suite
