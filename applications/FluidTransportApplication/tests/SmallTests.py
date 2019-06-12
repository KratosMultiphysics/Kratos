# Definition of the classes for the SMALL TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class steady_SC_test_2D(TF.TestFactory):
    file_name = "element_tests/steady_SC_test_2D/steady_SC_test_2D"
    file_parameters = "element_tests/steady_SC_test_2D/ProjectParameters.json"

class plume_test_2D(TF.TestFactory):
    file_name = "element_tests/plume_test_2D/pluma_2"
    file_parameters = "element_tests/plume_test_2D/ProjectParameters.json"

class conv_abso_diff_impl_NR_test_2D(TF.TestFactory):
    file_name = "element_tests/conv_abso_diff_impl_NR_test_2D/conv_abso_diff_impl_NR_test_2D"
    file_parameters = "element_tests/conv_abso_diff_impl_NR_test_2D/ProjectParameters.json"

class transient_wave_expl_cuad_test_2D(TF.TestFactory):
    file_name = "element_tests/transient_wave_expl_cuad_test_2D/transient_wave_expl_cuad_test_2D"
    file_parameters = "element_tests/transient_wave_expl_cuad_test_2D/ProjectParameters.json"

def SetTestSuite(suites):
    small_suite = suites['small']

    small_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            steady_SC_test_2D,
            plume_test_2D,
            conv_abso_diff_impl_NR_test_2D,
            transient_wave_expl_cuad_test_2D
        ])
    )

    return small_suite
