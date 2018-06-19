# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class Dam_Break_2D_Newtonian_Test(TF.TestFactory):
     file_name = "fluid_element_tests/Dam_Break_2D/Newtonian_fluid/Dam_break_2D"
     file_parameters = "fluid_element_tests/Dam_Break_2D/Newtonian_fluid/ProjectParameters.json"




def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            Dam_Break_2D_Newtonian_Test

        ])
    )

    return night_suite