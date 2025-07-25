# import Kratos

from KratosMultiphysics import *
from KratosMultiphysics.GeoMechanicsApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from test_excavation import KratosGeoMechanicsExcavationTests
from test_interface import KratosGeoMechanicsInterfaceTests
from test_reset_displacement import KratosGeoMechanicsResetDisplacementTests
from test_benchmark_set_1 import KratosGeoMechanicsBenchmarkSet1
from test_benchmark_set_2 import KratosGeoMechanicsBenchmarkSet2
from test_soil_structure_interactions import KratosGeoMechanicsSoilStructureInteractionTests
from test_water_pressure import KratosGeoMechanicsWaterPressureTests
from test_dynamics import KratosGeoMechanicsDynamicsTests
from test_dynamics_long import KratosGeoMechanicsDynamicsLongTests
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
from test_transient_thermal_validation import KratosGeoMechanicsTransientThermalValidationTests
from test_rotation_with_moving_load import KratosGeoMechanicsRotationWithMovingLoadTests
from test_time_integration import KratosGeoMechanicsTimeIntegrationTests
from c_phi_reduction_process import KratosGeoMechanicsCPhiReductionProcess
from test_partial_saturation import KratosGeoMechanicsPartialSaturation

from test_conditions import KratosGeoMechanicsConditionTests
from test_prescribed_derivatives import KratosGeoMechanicsPrescribedDerivatives
from test_dirichlet_u import KratosGeoMechanicsDirichletUTests
from test_normal_load_on_hexa_element import KratosGeoMechanicsNormalLoadHexaTests
from test_pressure_line_element import KratosGeoMechanicsTransientPressureLineElementTests
from test_pressure_point_flux import KratosGeoMechanicsTransientPressurePointFluxTests
from settlement_workflow import KratosGeoMechanicsSettlementWorkflowCppRoute, KratosGeoMechanicsSettlementWorkflowPyRoute
from test_compressibility import KratosGeoMechanicsCompressibilityTests
from fixed_spatial_variation import KratosGeoMechanicsFixedSpatialVariationTests
from test_integration_node_extrapolation import KratosGeoMechanicsExtrapolationTests
from test_truss_backbone_mat import KratosGeoMechanicsTrussBackboneMaterialTests
from test_line_interface_elements import KratosGeoMechanicsInterfaceElementTests
from test_three_dimensional_piping_validation import KratosGeoMechanicsThreeDimensionalPipingValidation
from test_master_slave_constraints import KratosGeoMechanicsMasterSlaveConstraints
from test_deactivation_with_structural_element import KratosGeoMechanicsDeactivationWithStructuralTest
from test_mohr_coulomb_with_tension_cutoff import KratosGeoMechanicsMohrCoulombWithTensionTests
from test_single_element_with_Mohr_Coulomb import KratosGeoMechanicsSingleElementWithMohrCoulomb
from one_dimensional_consolidation import KratosGeoMechanics1DConsolidation, KratosGeoMechanics1DConsolidationCppRoute
from test_apply_initial_uniform_stress_field import KratosGeoMechanicsApplyInitialUniformStressFieldTests
from test_dirichlet_release import KratosGeoMechanicsDirichletReleaseTests
from test_nodal_hydraulic_head import KratosGeoMechanicsHydraulicHeads
from test_set_moving_load_process import KratosGeoMechanicsSetMovingLoadProcessTests
from moving_column_with_fixed_pressure_above_phreatic_line import KratosGeoMechanicsMovingColumnWithFixedPressureAbovePhreaticLine

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
    small_test_cases = [
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
                        KratosGeoMechanicsConditionTests,
                        KratosGeoMechanicsPrescribedDerivatives,
                        KratosGeoMechanicsDirichletUTests,
                        KratosGeoMechanicsNormalLoadHexaTests,
                        KratosGeoMechanicsCompressibilityTests,
                        KratosGeoMechanicsFixedSpatialVariationTests,
                        KratosGeoMechanicsExtrapolationTests,
                        KratosGeoMechanicsTrussBackboneMaterialTests,
                        KratosGeoMechanicsInterfaceElementTests,
                        KratosGeoMechanicsMasterSlaveConstraints,
                        KratosGeoMechanicsSingleElementWithMohrCoulomb,
                        KratosGeoMechanicsMohrCoulombWithTensionTests,
                        KratosGeoMechanicsApplyInitialUniformStressFieldTests,
                        KratosGeoMechanicsDirichletReleaseTests,
                        KratosGeoMechanicsDeactivationWithStructuralTest,
                        KratosGeoMechanicsHydraulicHeads,
                        KratosGeoMechanicsSetMovingLoadProcessTests,
    ]

    night_test_cases = [
                        KratosGeoMechanicsSettlementWorkflowCppRoute,
                        KratosGeoMechanicsSettlementWorkflowPyRoute,
                        KratosGeoMechanicsCPhiReductionProcess,
                        KratosGeoMechanicsInterfaceTests,
                        KratosGeoMechanicsDynamicsTests,
                        KratosGeoMechanicsAbsorbingBoundaryColumnTests,
                        TestSellmeijersRule,
                        TestElementaryGroundWaterFlow,
                        KratosGeoMechanicsTransientThermalTests,
                        KratosGeoMechanicsTimeIntegrationTests,
                        KratosGeoMechanicsTransientPressureLineElementTests,
                        KratosGeoMechanicsTransientPressurePointFluxTests,
                        KratosGeoMechanicsPartialSaturation,
                        KratosGeoMechanics1DConsolidation,
                        KratosGeoMechanics1DConsolidationCppRoute,
                        KratosGeoMechanicsMovingColumnWithFixedPressureAbovePhreaticLine,
                        ]
    night_test_cases.extend(small_test_cases)

    # Create an array with all long tests only for validations
    valid_test_cases = [
                        KratosGeoMechanicsAbsorbingBoundaryColumn3DTests,
                        TestConsecutivePipeLines,
                        KratosGeoMechanicsBenchmarkSet1,
                        KratosGeoMechanicsBenchmarkSet2,
                        KratosGeoMechanicsTransientGroundWaterFlowTests,
                        TestSellmeijersRuleValidation,
                        KratosGeoMechanicsDynamicsLongTests,
                        KratosGeoMechanicsThreeDimensionalPipingValidation,
                        KratosGeoMechanicsTransientThermalValidationTests
                        ]

    # Create an array that contains all the tests from every testCase
    # in the list:

    all_test_cases = []
    all_test_cases.extend(night_test_cases)
    all_test_cases.extend(small_test_cases)
    all_test_cases.extend(valid_test_cases)
    suites = KratosUnittest.KratosSuites

    # add the tests to the corresponding suite,
    small_suite = suites['small']
    night_suite = suites['nightly']
    valid_suite = suites['validation']
    all_suite = suites['all']

    small_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(small_test_cases))
    night_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(night_test_cases))
    valid_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(valid_test_cases))
    all_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(all_test_cases))

    return suites


if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
