import subprocess
import os.path

# import Kratos
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.RANSApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utilities

# Import the tests o test_classes to create the suites
import run_cpp_unit_tests
from incompressible_potential_flow_solver_formulation_tests import IncompressiblePotentialFlowSolverFormulationTest

from monolithic_velocity_pressure_formulation_tests import MonolithicVelocityPressureFormulationTest
from monolithic_k_epsilon_high_re_formulation_tests import MonolithicKEpsilonHighReTest
from monolithic_k_omega_formulation_tests import MonolithicKOmegaTest
from monolithic_k_omega_sst_formulation_tests import MonolithicKOmegaSSTTest

from fractional_step_velocity_pressure_formulation_tests import FractionalStepVelocityPressureFormulationTest
from fractional_step_k_epsilon_high_re_formulation_tests import FractionalStepKEpsilonHighReTest
from fractional_step_k_omega_formulation_tests import FractionalStepKOmegaTest
from fractional_step_k_omega_sst_formulation_tests import FractionalStepKOmegaSSTTest

from custom_process_tests import CustomProcessTest

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([CustomProcessTest]))
    smallSuite.addTest(FractionalStepKOmegaSSTTest("testRfcVelocityTransient"))
    smallSuite.addTest(MonolithicKOmegaSSTTest("testRfcVelocityTransient"))

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # addign incompressible potential flow solver tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([IncompressiblePotentialFlowSolverFormulationTest]))

    # adding monolithic flow solver tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicVelocityPressureFormulationTest]))

    # adding monolithic k-epsilon high re tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKEpsilonHighReTest]))

    # adding monolithic k-omega tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKOmegaTest]))

    # adding monolithic k-omega-sst tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKOmegaSSTTest]))

    # adding fractional step flow solver tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepVelocityPressureFormulationTest]))

    # adding fractional step k-epsilon high re tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepKEpsilonHighReTest]))

    # adding fractional step k-omega tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepKOmegaTest]))

    # adding fractional step k-omega-sst tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepKOmegaSSTTest]))

    # For very long tests that should not be in nighly and you can use to validate
    # validationSuite = suites['validation']

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    # allSuite.addTests(validationSuite)

    return suites


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(
        KratosMultiphysics.Logger.Severity.WARNING)
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
