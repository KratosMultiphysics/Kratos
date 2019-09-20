# Definition of the classes for the SMALL TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

# class steady_SC_test_2D(TF.TestFactory):
#     file_name = "element_tests/steady_SC_test_2D/steady_SC_test_2D"
#     file_parameters = "element_tests/steady_SC_test_2D/ProjectParameters.json"

# class plume_test_2D(TF.TestFactory):
#     file_name = "element_tests/plume_test_2D/pluma_2"
#     file_parameters = "element_tests/plume_test_2D/ProjectParameters.json"

# class conv_abso_diff_impl_NR_test_2D(TF.TestFactory):
#     file_name = "element_tests/conv_abso_diff_impl_NR_test_2D/conv_abso_diff_impl_NR_test_2D"
#     file_parameters = "element_tests/conv_abso_diff_impl_NR_test_2D/ProjectParameters.json"

# class transient_wave_expl_cuad_test_2D(TF.TestFactory):
#     file_name = "element_tests/transient_wave_expl_cuad_test_2D/transient_wave_expl_cuad_test_2D"
#     file_parameters = "element_tests/transient_wave_expl_cuad_test_2D/ProjectParameters.json"

class SC_2D_test(TF.TestFactory):
    file_name = "element_tests/SC_2D_test/SC_2D_test"
    file_parameters = "element_tests/SC_2D_test/SC_2D_test.json"

class Q_source_2D_test(TF.TestFactory):
    file_name = "element_tests/Q_source_2D_test/Q_source_2D_test"
    file_parameters = "element_tests/Q_source_2D_test/Q_source_2D_test.json"

class Rot_2D_test(TF.TestFactory):
    file_name = "element_tests/Rot_2D_test/Rot_2D_test"
    file_parameters = "element_tests/Rot_2D_test/Rot_2D_test.json"

class DAC_2D_test(TF.TestFactory):
    file_name = "element_tests/DAC_2D_test/DAC_2D_test"
    file_parameters = "element_tests/DAC_2D_test/DAC_2D_test.json"

class Plume_PFEM2_2D_test(TF.TestFactory):
    file_name = "element_tests/Plume_PFEM2_2D_test/Plume_PFEM2_2D_test"
    file_parameters = "element_tests/Plume_PFEM2_2D_test/Plume_PFEM2_2D_test.json"

def SetTestSuite(suites):
    small_suite = suites['small']

    small_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            # steady_SC_test_2D,
            # plume_test_2D,
            # conv_abso_diff_impl_NR_test_2D,
            # transient_wave_expl_cuad_test_2D
            SC_2D_test,
            Q_source_2D_test,
            Rot_2D_test,
            DAC_2D_test,
            Plume_PFEM2_2D_test
        ])
    )

    return small_suite
