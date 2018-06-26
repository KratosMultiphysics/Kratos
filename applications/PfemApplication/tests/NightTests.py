# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class newtonian_dam_break_2D_test(TF.TestFactory):
     file_name = "fluid_tests/newtonian/dam_break_2D"
     file_parameters = None


def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
             newtonian_dam_break_2D_test
        ])
    )

    return night_suite
