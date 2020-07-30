# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest


class Bingham_Dam_Break_2D_Test(TF.TestFactory):
    file_name = "fluid_element_tests/Test_2D_Bingham/Test_2D_Bingham"
    file_parameters = "fluid_element_tests/Test_2D_Bingham/ProjectParameters.json"

class Water_Sloshing_3D_Test(TF.TestFactory):
    file_name = "fluid_element_tests/Test_3D_Newtonian_Sloshing/Test_3D_Newtonian_Sloshing"
    file_parameters = "fluid_element_tests/Test_3D_Newtonian_Sloshing/ProjectParameters.json"

class FSI_2D_Test(TF.TestFactory):
    file_name = "fluid_element_tests/Test_2D_FSI/Test_2D_FSI"
    file_parameters = "fluid_element_tests/Test_2D_FSI/ProjectParameters.json"


def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            Water_Sloshing_3D_Test,
        ])
    )

    return night_suite
