# Definition of the classes for the VALIDATION TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class Shear_Test_Von_Misses_Model(TF.TestFactory):
    file_materials  = "validation/shear_materials.json"
    file_parameters = "validation/shear_parameters.json"

class Shear_Traction_Test_Von_Misses_Model(TF.TestFactory):
    file_materials  = "validation/shear_materials.json"
    file_parameters = "validation/shear_traction_parameters.json"

class Shear_Test_Ogden_Model(TF.TestFactory):
    file_materials  = "validation/ogden_materials.json"
    file_parameters = "validation/shear_traction_parameters.json"

class Shear_Test_NeoHookean_Model(TF.TestFactory):
    file_materials  = "validation/neohookean_materials.json"
    file_parameters = "validation/shear_traction_parameters.json"

def SetTestSuite(suites):
    validation_suite = suites['validation']

    validation_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            Shear_Test_Von_Misses_Model,
            Shear_Traction_Test_Von_Misses_Model,
            Shear_Test_Ogden_Model,
            Shear_Test_NeoHookean_Model,
         ])
    )

    return validation_suite
