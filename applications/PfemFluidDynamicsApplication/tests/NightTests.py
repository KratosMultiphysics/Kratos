# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class Water_Sloshing_2D_Test(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_2D4N_shear"
    file_parameters = "element_tests/shear_2D_inputs.json"

def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            Water_Sloshing_2D_Test
        ])
    )

    return night_suite