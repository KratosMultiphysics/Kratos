# import Kratos

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from generalTests import KratosDamGeneralTests
from test_apply_load_vector_dam_processes import TestApplyLoadVectorDamProcesses
from test_apply_constraint_vector_dam_processes import (
    TestApplyConstraintVectorDamProcesses,
)
from test_impose_water_loads_condition_process import (
    TestImposeWaterLoadsConditionProcess,
)
from test_dam_linear_solver import TestDamLinearSolver
from test_dam_nodal_young_modulus_process import TestImposeNodalYoungModulusProcess
from test_dam_nodal_field_processes import TestDamNodalFieldProcesses
from test_impose_thermal_parameters_scalar_value_process import (
    TestImposeThermalParametersScalarValueProcess,
)
from test_dam_nonlocal_ownership import DamNonlocalOwnershipTest
from test_dam_process_lifecycle import TestDamProcessLifetime
from test_dam_analysis_process_lifecycle import TestDamAnalysisProcessLifecycle
from test_global_tangent_consistency import TestGlobalTangentConsistency
from test_interface_tangent_consistency import TestSmallDisplacementInterfaceTangent
from test_initial_gp_output_before_solve import TestInitialGaussPointOutputBeforeSolve
from test_dam_thermo_mechanic_initial_output import TestDamThermoMechanicInitialOutput
from test_dam_initial_stress_handling import TestDamInitialStressHandling
from test_dam_elastic_initial_state_transfer import TestDamElasticInitialStateTransfer
from test_dam_selective_initialization import TestDamSelectiveInitialization
from test_process_based_nodal_smoothing import DamProcessBasedNodalSmoothingTest
from test_dam_output_processes import TestDamOutputProcesses
from test_dam_gid_output_process import TestDamGiDOutputProcess
from test_dam_mapping_variables_utility import TestDamMappingVariablesUtility


def AssembleTestSuites():
    """Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    """

    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    # smallSuite will contain the following tests:
    # - testSmallExample
    smallSuite = suites["small"]
    smallSuite.addTest(KratosDamGeneralTests("test_joint_elastic_cohesive_2d_normal"))
    smallSuite.addTest(KratosDamGeneralTests("test_joint_elastic_cohesive_2d_shear"))
    smallSuite.addTest(
        KratosDamGeneralTests("test_joint_isotropic_damage_cohesive_2d_normal")
    )
    smallSuite.addTest(
        KratosDamGeneralTests("test_joint_isotropic_damage_cohesive_2d_shear")
    )
    smallSuite.addTest(KratosDamGeneralTests("test_joint_elastic_cohesive_3d_normal"))
    smallSuite.addTest(KratosDamGeneralTests("test_joint_elastic_cohesive_3d_shear"))
    smallSuite.addTest(
        KratosDamGeneralTests("test_joint_isotropic_damage_cohesive_3d_normal")
    )
    smallSuite.addTest(
        KratosDamGeneralTests("test_joint_isotropic_damage_cohesive_3d_shear")
    )
    smallSuite.addTest(KratosDamGeneralTests("test_construction"))
    smallSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases(
            [
                TestApplyLoadVectorDamProcesses,
                TestApplyConstraintVectorDamProcesses,
                TestImposeWaterLoadsConditionProcess,
                TestDamLinearSolver,
                TestGlobalTangentConsistency,
                TestImposeNodalYoungModulusProcess,
                TestDamNodalFieldProcesses,
                TestImposeThermalParametersScalarValueProcess,
                TestSmallDisplacementInterfaceTangent,
            ]
        )
    )
    smallSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases(
            [
                TestInitialGaussPointOutputBeforeSolve,
                TestDamThermoMechanicInitialOutput,
                TestDamInitialStressHandling,
                TestDamElasticInitialStateTransfer,
                DamProcessBasedNodalSmoothingTest,
                DamNonlocalOwnershipTest,
                TestDamProcessLifetime,
                TestDamAnalysisProcessLifecycle,
                TestDamOutputProcesses,
                TestDamGiDOutputProcess,
                TestDamMappingVariablesUtility,
            ]
        )
    )

    smallSuite.addTest(TestDamSelectiveInitialization("test_selective_gravity_stresses"))
    smallSuite.addTest(TestDamSelectiveInitialization("test_selective_gravity_and_temperature_increment"))
    smallSuite.addTest(TestDamSelectiveInitialization("test_rejects_unsupported_states_and_requires_gravity"))
    smallSuite.addTest(TestDamSelectiveInitialization("test_rejects_combination_with_previous_results_transfer"))
    smallSuite.addTest(TestDamSelectiveInitialization("test_does_not_modify_source_state"))
    smallSuite.addTest(TestDamSelectiveInitialization("test_destination_keeps_its_own_state_fields"))
    smallSuite.addTest(TestDamSelectiveInitialization("test_existing_reference_temperature_process"))
    smallSuite.addTest(TestDamSelectiveInitialization("test_use_process_based_local_equivalent_strain_is_registered"))

    # Create a test suit with the selected tests
    # nightSuite will contain the following tests:
    # - testSmallExample
    # - testNightlyFirstExample
    # - testNightlySecondExample
    nightSuite = suites["nightly"]
    nightSuite.addTests(smallSuite)

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites["all"]
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases(
            [
                KratosDamGeneralTests,
                TestApplyLoadVectorDamProcesses,
                TestApplyConstraintVectorDamProcesses,
                TestImposeWaterLoadsConditionProcess,
                TestDamLinearSolver,
                TestGlobalTangentConsistency,
                TestImposeNodalYoungModulusProcess,
                TestDamNodalFieldProcesses,
                TestImposeThermalParametersScalarValueProcess,
                TestInitialGaussPointOutputBeforeSolve,
                TestDamThermoMechanicInitialOutput,
                TestDamInitialStressHandling,
                TestDamElasticInitialStateTransfer,
                DamNonlocalOwnershipTest,
                DamProcessBasedNodalSmoothingTest,
                TestDamProcessLifetime,
                TestDamAnalysisProcessLifecycle,
                TestSmallDisplacementInterfaceTangent,
                TestDamOutputProcesses,
                TestDamGiDOutputProcess,
                TestDamMappingVariablesUtility,
            ]
        )
    )

    allSuite.addTest(TestDamSelectiveInitialization("test_selective_gravity_stresses"))
    allSuite.addTest(TestDamSelectiveInitialization("test_selective_gravity_and_temperature_increment"))
    allSuite.addTest(TestDamSelectiveInitialization("test_rejects_unsupported_states_and_requires_gravity"))
    allSuite.addTest(TestDamSelectiveInitialization("test_rejects_combination_with_previous_results_transfer"))
    allSuite.addTest(TestDamSelectiveInitialization("test_does_not_modify_source_state"))
    allSuite.addTest(TestDamSelectiveInitialization("test_destination_keeps_its_own_state_fields"))
    allSuite.addTest(TestDamSelectiveInitialization("test_existing_reference_temperature_process"))
    allSuite.addTest(TestDamSelectiveInitialization("test_use_process_based_local_equivalent_strain_is_registered"))

    return suites


if __name__ == "__main__":
    KratosUnittest.runTests(AssembleTestSuites())
