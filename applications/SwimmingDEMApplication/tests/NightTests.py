# Definition of the classes for the NIGHTLY TESTS

#Iimport Kratos
import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.SwimmingDEMApplication

# Import TestFactory
import InterpolationTestFactory as InterpolationTF
import TestFactory as TF
import FluidDEMTestFactory as FDEMTF
import SPFEMTestFactory as SPFEMTF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import importlib

class candelier_no_history_test(TF.TestFactory):
     file_name = "candelier_tests/candelier"
     file_parameters = "candelier_tests/ProjectParametersNoHistory.json"
     required_non_standard_modules = ['numpy', 'scipy']

class candelier_no_history_with_lift_test(TF.TestFactory):
     file_name = "candelier_tests/candelier"
     file_parameters = "candelier_tests/ProjectParametersNoHistoryWithLift.json"
     required_non_standard_modules = ['numpy', 'scipy']

class candelier_no_history_non_inertial_test(TF.TestFactory):
     file_name = "candelier_tests/candelier"
     file_parameters = "candelier_tests/ProjectParametersNoHistoryNonInertial.json"
     required_non_standard_modules = ['numpy', 'scipy']

class candelier_with_history_test(TF.TestFactory):
     file_name = "candelier_tests/candelier"
     file_parameters = "candelier_tests/ProjectParametersWithHistory.json"
     required_non_standard_modules = ['numpy', 'scipy']

class candelier_with_history_hinsberg_test(TF.TestFactory):
     file_name = "candelier_tests/candelier"
     file_parameters = "candelier_tests/ProjectParametersWithHistoryHinsberg.json"
     required_non_standard_modules = ['numpy', 'scipy']

# This test is ready to run but the implementation is not complete
# # (it is non-trivial), so the result is not correct
# class candelier_with_history_non_inertial_test(TF.TestFactory):
#      file_name = "candelier_tests/candelier"
#      file_parameters = "candelier_tests/ProjectParametersWithHistoryNonInertial.json"
#      required_non_standard_modules = ['numpy', 'scipy']

class interpolation_test_linear(InterpolationTF.TestFactory):
     file_name = "interpolation_tests/cube"
     file_parameters = "interpolation_tests/ProjectParametersCubeLinear.json"

class interpolation_test_nonlinear_time_no_substepping(InterpolationTF.TestFactory):
     file_name = "interpolation_tests/cube"
     file_parameters = "interpolation_tests/ProjectParametersCubeNonlinearTimeNoSubstepping.json"


# List of all classes above
list_of_classes = [test_class for test_class in
                    TF.TestFactory.__subclasses__()
                  + InterpolationTF.TestFactory.__subclasses__()]

def CheckImportsArePossible(test):
     if hasattr(test, 'required_non_standard_modules'):
          for module in test.required_non_standard_modules:
               try:
                    importlib.import_module(module)
               except:
                    return False
     return True

def SetTestSuite(suites):
    night_suite = suites['nightly']

    available_tests = [test for test in list_of_classes if CheckImportsArePossible(test)]

    night_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(available_tests))

    return night_suite

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites
    night_suite = SetTestSuite(suites)
    suites['all'].addTests(night_suite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.DETAIL)
    KratosUnittest.runTests(AssembleTestSuites())
