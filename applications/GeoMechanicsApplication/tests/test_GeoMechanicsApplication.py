# import Kratos

from KratosMultiphysics import *
from KratosMultiphysics.GeoMechanicsApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

import run_cpp_unit_tests

# Import the tests o test_classes to create the suits
from generalTests import KratosGeoMechanicsGeneralTests
from test_excavation import KratosGeoMechanicsExcavationTests
from test_interface import KratosGeoMechanicsInterfaceTests
from test_reset_displacement import KratosGeoMechanicsResetDisplacementTests
from test_benchmark_set_1 import KratosGeoMechanicsBenchmarkSet1
from test_benchmark_set_2 import KratosGeoMechanicsBenchmarkSet2
from test_soil_structure_interactions import KratosGeoMechanicsSoilStructureInteractionTests
from test_water_pressure import KratosGeoMechanicsWaterPressureTests
from test_dynamics import KratosGeoMechanicsDynamicsTests
from test_elements import KratosGeoMechanicsElementTypeTests
from test_steady_state_groundwater_flow import KratosGeoMechanicsSteadyStateGroundWaterFlowTests
from test_transient_groundwater_flow import KratosGeoMechanicsTransientGroundWaterFlowTests
from test_soil_weight import KratosGeoMechanicsSoilWeightTests
from test_settlement import KratosGeoMechanicsSettlementTests
from test_curved_beam_elements import KratosGeoMechanicsCurvedBeamElementTests
from absorbing_boundary import KratosGeoMechanicsAbsorbingBoundaryColumnTests
from absorbing_boundary_3D import KratosGeoMechanicsAbsorbingBoundaryColumn3DTests
from test_elementary_groundwater_flow import TestElementaryGroundWaterFlow
from test_sellmeijers_rule import TestSellmeijersRule
from test_sellmeijers_rule_validation import TestSellmeijersRuleValidation
from test_consecutive_pipe_lines import TestConsecutivePipeLines
from test_line_loads import KratosGeoMechanicsLineLoadTests
from test_element_lab import KratosGeoMechanicsLabElementTests
from test_parameter_field import KratosGeoMechanicsParameterFieldTests
from test_normal_load_on_1d_element import KratosGeoMechanicsNormalLoad1DTests
from test_k0_procedure_process import KratosGeoMechanicsK0ProcedureProcessTests
from test_geomechanics_solver import KratosGeoMechanicsSolverTests
from test_column_changing_waterlevel import KratosGeoMechanicsChangingWaterLevelTests
from test_set_multiple_moving_load_process import KratosGeoMechanicsSetMultipleMovingLoadProcessTests
from test_strain_measures import KratosGeoMechanicsStrainMeasureTests
from test_transient_thermal import KratosGeoMechanicsTransientThermalTests
from test_rotation_with_moving_load import KratosGeoMechanicsRotationWithMovingLoadTests
from test_time_integration import KratosGeoMechanicsTimeIntegrationTests
from test_conditions import KratosGeoMechanicsConditionTests

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    # Create an array with the selected tests (Small tests):
    # smallSuite will contain the following tests:
    # - testSmallExample

    small_test_cases = [
                        KratosGeoMechanicsGeneralTests,
                        KratosGeoMechanicsExcavationTests,
                        KratosGeoMechanicsResetDisplacementTests,
                        KratosGeoMechanicsSoilStructureInteractionTests,
                        KratosGeoMechanicsWaterPressureTests,
                        KratosGeoMechanicsElementTypeTests,
                        KratosGeoMechanicsSteadyStateGroundWaterFlowTests,
                        KratosGeoMechanicsSoilWeightTests,
                        KratosGeoMechanicsSettlementTests,
                        KratosGeoMechanicsLineLoadTests,
                        KratosGeoMechanicsCurvedBeamElementTests,
                        KratosGeoMechanicsLabElementTests,
                        KratosGeoMechanicsParameterFieldTests,
                        KratosGeoMechanicsNormalLoad1DTests,
                        KratosGeoMechanicsK0ProcedureProcessTests,
                        KratosGeoMechanicsSolverTests,
                        KratosGeoMechanicsChangingWaterLevelTests,
                        KratosGeoMechanicsStrainMeasureTests,
                        KratosGeoMechanicsSetMultipleMovingLoadProcessTests,
                        KratosGeoMechanicsRotationWithMovingLoadTests,
                        KratosGeoMechanicsConditionTests
                        ]

    # Create an array with the selected tests
    # nightSuite will contain the following tests:
    # - testSmallExample
    # - testNightlyFirstExample
    # - testNightlySecondExample

    night_test_cases = [
                        KratosGeoMechanicsInterfaceTests,
                        KratosGeoMechanicsDynamicsTests,
                        KratosGeoMechanicsAbsorbingBoundaryColumnTests,
                        TestSellmeijersRule,
                        TestElementaryGroundWaterFlow,
                        KratosGeoMechanicsTransientThermalTests,
                        KratosGeoMechanicsTimeIntegrationTests
                        ]
    night_test_cases.extend(small_test_cases)

    # Create an array with all long tests only for validations
    valid_test_cases = [
                        KratosGeoMechanicsAbsorbingBoundaryColumn3DTests,
                        TestConsecutivePipeLines,
                        KratosGeoMechanicsBenchmarkSet1,
                        KratosGeoMechanicsBenchmarkSet2,
                        KratosGeoMechanicsTransientGroundWaterFlowTests,
                        TestSellmeijersRuleValidation
                        ]

    # Create an array that contains all the tests from every testCase
    # in the list:

    all_test_cases = []
    all_test_cases.extend(night_test_cases)
    all_test_cases.extend(small_test_cases)
    all_test_cases.extend(valid_test_cases)
    suites = KratosUnittest.KratosSuites

    # add the tests to the corresponding suite,
    smallSuite = suites['small']
    nightSuite = suites['nightly']
    validSuite = suites['validation']
    allSuite = suites['all']

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(small_test_cases))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(night_test_cases))
    validSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(valid_test_cases))
    allSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(all_test_cases))

    return suites


if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")
