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

# process test_classes
from custom_process_tests import CustomProcessTest

# flow solver test_classes
from incompressible_potential_flow_solver_formulation_tests import IncompressiblePotentialFlowSolverFormulationTest
from monolithic_velocity_pressure_formulation_tests import MonolithicVelocityPressureFormulationTest
from fractional_step_velocity_pressure_formulation_tests import FractionalStepVelocityPressureFormulationTest

# turbulence model test_classes
### k-epsilon test_classes
from monolithic_k_epsilon_formulation_tests import MonolithicKEpsilonTest
from monolithic_k_epsilon_formulation_tests import MonolithicKEpsilonPeriodicTest
from fractional_step_k_epsilon_formulation_tests import FractionalStepKEpsilonTest

### k-omega test_classes
from monolithic_k_omega_formulation_tests import MonolithicKOmegaTest
from monolithic_k_omega_formulation_tests import MonolithicKOmegaPeriodicTest
from fractional_step_k_omega_formulation_tests import FractionalStepKOmegaTest

### k-omega-sst test_classes
from monolithic_k_omega_sst_formulation_tests import MonolithicKOmegaSSTTest
from monolithic_k_omega_sst_formulation_tests import MonolithicKOmegaSSTPeriodicTest
from fractional_step_k_omega_sst_formulation_tests import FractionalStepKOmegaSSTTest

### adjoint two element test_classes
from adjoint_cc_two_elements_tests import AdjointCircularConvectionTwoElementsTest
from adjoint_diffusion_two_elements_tests import AdjointDiffusionTwoElementsTest
from adjoint_ke_two_elements_tests import AdjointKEpsilonTwoElementsTest
from adjoint_kw_two_elements_tests import AdjointKOmegaTwoElementsTest

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

    # adding process tests to small suite
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([CustomProcessTest]))

    # adding k-epsilon high re periodic tests to small suite
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKEpsilonPeriodicTest]))

    # adding k-omega high re periodic tests to small suite
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKOmegaPeriodicTest]))

    # adding k-omega-sst high re periodic tests to small suite
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKOmegaSSTPeriodicTest]))

    # adding representative transient tests to small suite
    smallSuite.addTest(FractionalStepKOmegaSSTTest("testVMSRfcVelocityTransient"))
    smallSuite.addTest(MonolithicKOmegaSSTTest("testVMSRfcVelocityTransient"))
    smallSuite.addTest(MonolithicKOmegaSSTTest("testQSVMSRfcVelocityTransient"))

    # adding adjoint two element tests
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([AdjointKEpsilonTwoElementsTest]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([AdjointCircularConvectionTwoElementsTest]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([AdjointDiffusionTwoElementsTest]))

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # adding incompressible potential flow solver tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([IncompressiblePotentialFlowSolverFormulationTest]))

    # adding monolithic flow solver tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicVelocityPressureFormulationTest]))

    # adding fractional step flow solver tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepVelocityPressureFormulationTest]))

    # adding monolithic k-epsilon high re tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKEpsilonTest]))
    # adding fractional step k-epsilon high re tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepKEpsilonTest]))

    # adding monolithic k-omega high re tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicKOmegaTest]))
    # adding fractional step k-omega high re tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepKOmegaTest]))

    # adding monolithic k-omega-sst high re tests
    nightSuite.addTest(MonolithicKOmegaSSTTest("testVMSAfcTkeSteady"))
    nightSuite.addTest(MonolithicKOmegaSSTTest("testVMSAfcVelocitySteady"))
    nightSuite.addTest(MonolithicKOmegaSSTTest("testVMSRfcTkeSteady"))
    nightSuite.addTest(MonolithicKOmegaSSTTest("testVMSRfcVelocitySteady"))
    nightSuite.addTest(MonolithicKOmegaSSTTest("testVMSRfcTkeTransient"))
    nightSuite.addTest(MonolithicKOmegaSSTTest("testQSVMSRfcVelocitySteady"))

    # adding fractional step k-omega-sst high re tests
    nightSuite.addTest(FractionalStepKOmegaSSTTest("testVMSAfcTkeSteady"))
    nightSuite.addTest(FractionalStepKOmegaSSTTest("testVMSAfcVelocitySteady"))
    nightSuite.addTest(FractionalStepKOmegaSSTTest("testVMSRfcTkeSteady"))
    nightSuite.addTest(FractionalStepKOmegaSSTTest("testVMSRfcVelocitySteady"))
    nightSuite.addTest(FractionalStepKOmegaSSTTest("testVMSRfcTkeTransient"))

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([AdjointKOmegaTwoElementsTest]))

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)
    allSuite.addTests(validationSuite)

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

