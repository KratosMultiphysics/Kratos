# Definition of the classes for the SMALL TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class consolidation_2D(TF.TestFactory):
    file_name = "element_tests/consolidation_2D/consolidation_2D"
    file_parameters = "element_tests/consolidation_2D/ProjectParameters.json"

class consolidation_interface_2D(TF.TestFactory):
    file_name = "element_tests/consolidation_interface_2D/consolidation_interface_2D"
    file_parameters = "element_tests/consolidation_interface_2D/ProjectParameters.json"

class interface_elastic_linear(TF.TestFactory):
    file_name = "constitutive_law_tests/interface_elastic_linear/interface_elastic_linear"
    file_parameters = "constitutive_law_tests/interface_elastic_linear/ProjectParameters.json"

class interface_isotropic_damage(TF.TestFactory):
    file_name = "constitutive_law_tests/interface_isotropic_damage_linear/Test_UnixialTension"
    file_parameters = "constitutive_law_tests/interface_isotropic_damage_linear/ProjectParameters.json"

class interface_mc_tension_cutoff(TF.TestFactory):
    file_name = "constitutive_law_tests/interface_MC_tension_cutoff/Test_MixedMode3D"
    file_parameters = "constitutive_law_tests/interface_MC_tension_cutoff//ProjectParameters.json"

class nodal_discharge_2D(TF.TestFactory):
    file_name = "element_tests/nodal_discharge/2_NodalDischarge3elem"
    file_parameters = "element_tests/nodal_discharge/ProjectParameters.json"

def SetTestSuite(suites):
    small_suite = suites['small']

    small_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            consolidation_2D,
            consolidation_interface_2D,
            interface_elastic_linear,
            interface_isotropic_damage,
            interface_mc_tension_cutoff,
            nodal_discharge_2D
        ])
    )

    return small_suite
