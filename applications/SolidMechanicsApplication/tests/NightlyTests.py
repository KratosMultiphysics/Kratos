# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF
    
#class ValidationTest(TestFactory):
    #file_name = "path_to_my_test"
    
class SmallStrains_IsotropicDamage_SimoJu_Test(TF.TestFactory):
    file_name = "material_models_tests/isotropic_damage_material_model/PlaneStress_FourPointShear_test"
