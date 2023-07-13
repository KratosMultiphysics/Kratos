# Definition of the classes for the SMALL TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class undrained_soil_column_2D(TF.TestFactory):
    file_name = "element_tests/undrained_soil_column_2D/undrained_soil_column_2D"
    file_parameters = "element_tests/undrained_soil_column_2D/ProjectParameters.json"

class interface_elastic_linear(TF.TestFactory):
    file_name = "constitutiveModel_tests/interface_elastic_linear/interface_elastic_linear"
    file_parameters = "constitutiveModel_tests/interface_elastic_linear/ProjectParameters.json"

class interface_isotropic_damage(TF.TestFactory):
    file_name = "constitutiveModel_tests/interface_isotropic_damage_linear/Test_UnixialTension"
    file_parameters = "constitutiveModel_tests/interface_isotropic_damage_linear/ProjectParameters.json"

def SetTestSuite(suites):
    small_suite = suites['small']

    small_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            undrained_soil_column_2D,
            interface_elastic_linear,
            interface_isotropic_damage
        ])
    )

    return small_suite
