# Definition of the classes for the SMALL TESTS

#Iimport Kratos
import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.SwimmingDEMApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Importing test factories if possible

try:
     import InterpolationTestFactory as InterpolationTF
     interpolation_imports_available = True
except ImportError:
     interpolation_imports_available = False
try:
     import CandelierTestFactory as CandelierTF
     candelier_imports_available = True
except ImportError:
     candelier_imports_available = False
try:
     import FluidDEMTestFactory as FDEMTF
     fluid_DEM_coupling_imports_available = True
except ImportError:
     fluid_DEM_coupling_imports_available = False
# List of tests that are available
available_tests = []

if interpolation_imports_available:
     class interpolation_test_linear(InterpolationTF.TestFactory):
          file_name = "interpolation_tests/cube"
          file_parameters = "interpolation_tests/ProjectParametersCubeLinear.json"

     class interpolation_test_nonlinear_time_no_substepping(InterpolationTF.TestFactory):
          file_name = "interpolation_tests/cube"
          file_parameters = "interpolation_tests/ProjectParametersCubeNonlinearTimeNoSubstepping.json"

     available_tests += [test_class for test_class in InterpolationTF.TestFactory.__subclasses__()]


if candelier_imports_available:
     class candelier_no_history_test(CandelierTF.TestFactory):
          file_name = "candelier_tests/candelier"
          file_parameters = "candelier_tests/ProjectParametersNoHistory.json"

     class candelier_no_history_with_lift_test(CandelierTF.TestFactory):
          file_name = "candelier_tests/candelier"
          file_parameters = "candelier_tests/ProjectParametersNoHistoryWithLift.json"

     class candelier_no_history_non_inertial_test(CandelierTF.TestFactory):
          file_name = "candelier_tests/candelier"
          file_parameters = "candelier_tests/ProjectParametersNoHistoryNonInertial.json"

     class candelier_with_history_test(CandelierTF.TestFactory):
          file_name = "candelier_tests/candelier"
          file_parameters = "candelier_tests/ProjectParametersWithHistory.json"

     class candelier_with_history_hinsberg_test(CandelierTF.TestFactory):
          file_name = "candelier_tests/candelier"
          file_parameters = "candelier_tests/ProjectParametersWithHistoryHinsberg.json"

     # This test is ready to run but the implementation is not complete
     # # (it is non-trivial), so the result is not correct
     # class candelier_with_history_non_inertial_test(CandelierTF.TestFactory):
     #      file_name = "candelier_tests/candelier"
     #      file_parameters = "candelier_tests/ProjectParametersWithHistoryNonInertial.json"

     available_tests += [test_class for test_class in CandelierTF.TestFactory.__subclasses__()]

if fluid_DEM_coupling_imports_available:
     class fluid_dem_coupling_one_way_test(FDEMTF.TestFactory):
          file_name = "fluid_dem_tests/settling_cube"
          file_parameters = "fluid_dem_tests/ProjectParameters.json"

     available_tests += [test_class for test_class in FDEMTF.TestFactory.__subclasses__()]

def SetTestSuite(suites):
    small_suite = suites['small']

    small_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(available_tests))

    return small_suite

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites
    small_suite = SetTestSuite(suites)
    suites['all'].addTests(small_suite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())

