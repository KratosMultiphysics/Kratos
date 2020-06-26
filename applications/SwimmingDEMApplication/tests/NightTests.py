# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import InterpolationTestFactory as InterpolationTF
import TestFactory as TF
import FluidDEMTestFactory as FDEMTF
import SPFEMTestFactory as SPFEMTF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class candelier_no_history_test(TF.TestFactory):
     file_name = "candelier_tests/candelier"
     file_parameters = "candelier_tests/ProjectParametersNoHistory.json"

class candelier_no_history_with_lift_test(TF.TestFactory):
     file_name = "candelier_tests/candelier"
     file_parameters = "candelier_tests/ProjectParametersNoHistoryWithLift.json"

class candelier_no_history_non_inertial_test(TF.TestFactory):
     file_name = "candelier_tests/candelier"
     file_parameters = "candelier_tests/ProjectParametersNoHistoryNonInertial.json"

class candelier_with_history_test(TF.TestFactory):
     file_name = "candelier_tests/candelier"
     file_parameters = "candelier_tests/ProjectParametersWithHistory.json"

class candelier_with_history_hinsberg_test(TF.TestFactory):
     file_name = "candelier_tests/candelier"
     file_parameters = "candelier_tests/ProjectParametersWithHistoryHinsberg.json"

# This test is ready to run but the implementation is not complete
# (it is non-trivial), so the result is not correct
class candelier_with_history_non_inertial_test(TF.TestFactory):
     file_name = "candelier_tests/candelier"
     file_parameters = "candelier_tests/ProjectParametersWithHistoryNonInertial.json"

class interpolation_test_linear(InterpolationTF.TestFactory):
     file_name = "interpolation_tests/cube"
     file_parameters = "interpolation_tests/ProjectParametersCubeLinear.json"

class interpolation_test_nonlinear_time_no_substepping(InterpolationTF.TestFactory):
     file_name = "interpolation_tests/cube"
     file_parameters = "interpolation_tests/ProjectParametersCubeNonlinearTimeNoSubstepping.json"

class fluid_dem_coupling_one_way_test(FDEMTF.TestFactory):
     file_name = "fluid_dem_tests/settling_cube"
     file_parameters = "fluid_dem_tests/ProjectParameters.json"

class sdem_pfem_coupling_one_way_test(SPFEMTF.TestFactory):
     file_name = "PFEM-DEM_tests/sdem_pfem_coupling_one_way_test"
     file_parameters = "PFEM-DEM_tests/ProjectParameters.json"

def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
          candelier_no_history_test,
          candelier_no_history_with_lift_test,
          candelier_no_history_non_inertial_test,
          candelier_with_history_test,
          candelier_with_history_hinsberg_test,
          interpolation_test_linear,
          interpolation_test_nonlinear_time_no_substepping,
          fluid_dem_coupling_one_way_test,
          sdem_pfem_coupling_one_way_test
          ])
    )

    return night_suite
