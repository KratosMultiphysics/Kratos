# Definition of the classes for the SMALL TESTS

# Import Kratos and necessary applications
import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.SwimmingDEMApplication
from KratosMultiphysics import Logger

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as UnitTest
import BackwardCouplingTestFactory as BackwardCouplingTF

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
     Logger.PrintWarning("SwimmingDEMTests", "Failed to import some of the modules necessary for the Candelier tests.")
try:
     import FluidDEMTestFactory as FDEMTF
     fluid_DEM_coupling_imports_available = True
except ImportError:
     fluid_DEM_coupling_imports_available = False
try:
     import AnalyticTestFactory as AnalyticTF
     analytic_imports_available = True
except ImportError:
     analytic_imports_available = False

class interpolation_test_linear(InterpolationTF.TestFactory):
     file_name = "interpolation_tests/cube"
     file_parameters = "interpolation_tests/ProjectParametersCubeLinear.json"

class interpolation_test_nonlinear_time_no_substepping(InterpolationTF.TestFactory):
     file_name = "interpolation_tests/cube"
     file_parameters = "interpolation_tests/ProjectParametersCubeNonlinearTimeNoSubstepping.json"

class backward_coupling_single_particle_no_time_filter(BackwardCouplingTF.TestFactory):
     file_name = "backward_coupling_tests/cube_single_particle"
     file_parameters = "backward_coupling_tests/ProjectParametersCubeNoTimeFilter.json"

class backward_coupling_two_balls_periodic_no_time_filter(BackwardCouplingTF.TestFactory):
     file_name = "backward_coupling_tests/cube_two_balls_periodic"
     file_parameters = "backward_coupling_tests/ProjectParametersCubeTwoBallsPeriodic.json"

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

class fluid_dem_coupling_one_way_test(FDEMTF.TestFactory):
     file_name = "fluid_dem_tests/settling_cube"
     file_parameters = "fluid_dem_tests/ProjectParameters.json"

class CFD_DEM_two_way_test(FDEMTF.TestFactory):
     file_name = "CFD_DEM_two_way_tests/Two_way_testFluid"
     file_parameters = "CFD_DEM_two_way_tests/ProjectParameters.json"

class CFD_DEM_two_way_test_DVMS(FDEMTF.TestFactory):
     file_name = "CFD_DEM_two_way_tests/Two_way_testFluid_DVMS"
     file_parameters = "CFD_DEM_two_way_tests/ProjectParameters_DVMS.json"

class analytic_multiple_ghosts_test(AnalyticTF.TestFactory):
     file_name = "analytic_tests/multiple_phantoms/multiple_phantom_test"
     file_parameters = "analytic_tests/multiple_phantoms/ProjectParameters.json"

available_tests = []
available_tests += [test_class for test_class in InterpolationTF.TestFactory.__subclasses__()]
available_tests += [test_class for test_class in BackwardCouplingTF.TestFactory.__subclasses__()]
if candelier_imports_available:
     available_tests += [test_class for test_class in CandelierTF.TestFactory.__subclasses__()]
available_tests += [test_class for test_class in FDEMTF.TestFactory.__subclasses__()]
available_tests += [test_class for test_class in AnalyticTF.TestFactory.__subclasses__()]

def SetTestSuite(suites):
    small_suite = suites['small']
    small_suite.addTests(UnitTest.TestLoader().loadTestsFromTestCases(available_tests))

    return small_suite

def AssembleTestSuites():
    suites = UnitTest.KratosSuites
    small_suite = SetTestSuite(suites)
    suites['all'].addTests(small_suite)

    return suites

if __name__ == '__main__':
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    UnitTest.runTests(AssembleTestSuites())
    UnitTest.main()

