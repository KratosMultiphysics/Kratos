# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF
    

class SmallStrains_IsotropicDamage_SimoJu_Test(TF.TestFactory):
    file_name = "material_tests/isotropic_damage_material_model/PlaneStress_FourPointShear"
    file_parameters = None
