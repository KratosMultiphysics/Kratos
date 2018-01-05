# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF
    
# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class SmallStrains_IsotropicDamage_SimoJu_Test(TF.TestFactory):
    file_name = "material_tests/isotropic_damage_material_model/PlaneStress_FourPointShear"
    file_parameters = None

def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            SmallStrains_IsotropicDamage_SimoJu_Test
        ])
    )
    
    return night_suite
