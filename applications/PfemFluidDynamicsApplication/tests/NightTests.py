# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class Dam_Break_2D_Newtonian_Test(TF.TestFactory):
     file_name = "fluid_element_tests/Dam_Break_2D/Newtonian_fluid/Dam_break_2D"
     file_parameters = "fluid_element_tests/Dam_Break_2D/Newtonian_fluid/ProjectParameters.json"
class Water_sloshing_Box_3D_Newtonian_Test(TF.TestFactory):
     file_name = "fluid_element_tests/Water_sloshing_Box/Newtonian_fluid/Water_sloshing_Box"
     file_parameters = "fluid_element_tests/Water_sloshing_Box/Newtonian_fluid/ProjectParameters.json"
class Water_sloshing_Box_3D_Non_Newtonian_Test(TF.TestFactory):
      file_name = "fluid_element_tests/Water_sloshing_Box/Non_Newtonian_fluid/Water_sloshing_Box"
      file_parameters = "fluid_element_tests/Water_sloshing_Box/Non_Newtonian_fluid/ProjectParameters.json"
class Bingham_Dam_Break_2D_Test(TF.TestFactory):
    file_name = "fluid_element_tests/Test_2D_Bingham/Test_2D_Bingham"
    file_parameters = "fluid_element_tests/Test_2D_Bingham/ProjectParameters.json"



def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            # Dam_Break_2D_Newtonian_Test,
            # Water_sloshing_Box_3D_Non_Newtonian_Test,
            # Water_sloshing_Box_3D_Newtonian_Test
            Bingham_Dam_Break_2D_Test
        ])
    )

    return night_suite
