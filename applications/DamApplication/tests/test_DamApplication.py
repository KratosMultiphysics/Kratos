# import Kratos

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from generalTests import KratosDamGeneralTests
from test_apply_load_vector_dam_processes import TestApplyLoadVectorDamProcesses
from test_dam_linear_solver import TestDamLinearSolver
from test_dam_nodal_young_modulus_process import TestImposeNodalYoungModulusProcess
from test_dam_nonlocal_ownership import DamNonlocalOwnershipTest
from test_dam_process_lifecycle import TestDamProcessLifetime
from test_global_tangent_consistency import TestGlobalTangentConsistency
from test_initial_gp_output_before_solve import TestInitialGaussPointOutputBeforeSolve
from test_interface_tangent_consistency import TestSmallDisplacementInterfaceTangent
from test_process_based_nodal_smoothing import DamProcessBasedNodalSmoothingTest


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
                TestDamLinearSolver,
                TestGlobalTangentConsistency,
                TestImposeNodalYoungModulusProcess,
                TestSmallDisplacementInterfaceTangent,
            ]
        )
    )
    smallSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases(
            [
                TestApplyLoadVectorDamProcesses,
                TestInitialGaussPointOutputBeforeSolve,
                DamProcessBasedNodalSmoothingTest,
                DamNonlocalOwnershipTest,
                TestDamProcessLifetime,
            ]
        )
    )

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
                TestDamLinearSolver,
                TestGlobalTangentConsistency,
                TestImposeNodalYoungModulusProcess,
                TestInitialGaussPointOutputBeforeSolve,
                TestDamProcessLifetime,
                DamProcessBasedNodalSmoothingTest,
                DamNonlocalOwnershipTest,
                TestSmallDisplacementInterfaceTangent,
            ]
        )
    )

    return suites


if __name__ == "__main__":
    KratosUnittest.runTests(AssembleTestSuites())
