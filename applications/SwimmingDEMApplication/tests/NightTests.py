# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF
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

class Test(SPFEMTF.TestFactory):
     file_name = "PFEM-DEM_tests/Test"
     file_parameters = "PFEM-DEM_tests/ProjectParameters.json"

# This test is ready to run but the implementation is not complete
# (it is non-trivial), so the result is not correct
class candelier_with_history_non_inertial_test(TF.TestFactory):
     file_name = "candelier_tests/candelier"
     file_parameters = "candelier_tests/ProjectParametersWithHistoryNonInertial.json"

def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
          candelier_no_history_test,
          candelier_no_history_with_lift_test,
          candelier_no_history_non_inertial_test,
          candelier_with_history_test,
          candelier_with_history_hinsberg_test,
          Test
          ])
    )

    return night_suite
