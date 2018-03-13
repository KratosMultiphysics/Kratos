# Definition of the classes for the VALIDATION TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class Dynamic_Test_Bossak_TL_3D(TF.TestFactory):
    file_name = "dynamic_tests/solid_elements/dynamic_bossak_TL3D"
    file_parameters = None

class Dynamic_Test_Simo_TL_3D(TF.TestFactory):
    file_name = "dynamic_tests/solid_elements/dynamic_bossak_TL3D"
    file_parameters  = "dynamic_tests/solid_elements/dynamic_simo_TL3D_input.json"



def SetTestSuite(suites):
    validation_suite = suites['validation']

    validation_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            Dynamic_Test_Bossak_TL_3D,
            Dynamic_Test_Simo_TL_3D
        ])
    )

    return validation_suite
