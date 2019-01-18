# Definition of the classes for the SMALL TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class Water_Sloshing_2D_Test_Newtonian(TF.TestFactory):
    file_name = "fluid_element_tests/Water_sloshing_2D/Newtonian_fluid/Water_sloshing_2D"
    file_parameters = "fluid_element_tests/Water_sloshing_2D/Newtonian_fluid/ProjectParameters.json"

class Water_Sloshing_2D_Test_Non_Newtonian(TF.TestFactory):
    file_name = "fluid_element_tests/Water_sloshing_2D/Non_Newtonian_fluid/Water_sloshing_2D"
    file_parameters = "fluid_element_tests/Water_sloshing_2D/Non_Newtonian_fluid/ProjectParameters.json"

class Water_Sloshing_2D_Test_Mu_Rheology(TF.TestFactory):
    file_name = "fluid_element_tests/Water_sloshing_2D/Mu_Rhelogy/Water_sloshing_2D"
    file_parameters = "fluid_element_tests/Water_sloshing_2D/Mu_Rheology/ProjectParameters.json"

class MuIrheology_Dam_Break_2D_Test(TF.TestFactory):
    file_name = "fluid_element_tests/Test_2D_muIrheology/Test_2D_muIrheology"
    file_parameters = "fluid_element_tests/Test_2D_muIrheology/ProjectParameters.json"

class Water_Sloshing_2D_Test(TF.TestFactory):
    file_name = "fluid_element_tests/Test_2D_Newtonian_Sloshing/Test_2D_Newtonian_Sloshing"
    file_parameters = "fluid_element_tests/Test_2D_Newtonian_Sloshing/ProjectParameters.json"


def SetTestSuite(suites):
    small_suite = suites['small']

    small_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
           # Water_Sloshing_2D_Test_Newtonian,
           # Water_Sloshing_2D_Test_Non_Newtonian,
           # Water_Sloshing_2D_Test_Mu_Rheology,
           MuIrheology_Dam_Break_2D_Test,
           Water_Sloshing_2D_Test
        ])
    )

    return small_suite
