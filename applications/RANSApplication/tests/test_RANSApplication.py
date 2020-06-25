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

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    nightSuite.addTest((IncompressiblePotentialFlowSolverFormulationTest('testIncompressiblePotentialFlow')))
    nightSuite.addTest((MonolithicVelocityPressureFormulationTest('testMonolithicVelocityPressureSteady')))
    nightSuite.addTest((MonolithicVelocityPressureFormulationTest('testMonolithicVelocityPressureTransient')))

    # adding monolithic k-epsilon high re tests
    nightSuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReAfcTkeSteady')))
    nightSuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReAfcVelocitySteady')))
    nightSuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReRfcTkeSteady')))
    nightSuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReRfcVelocitySteady')))
    nightSuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReRfcTkeTransient')))
    nightSuite.addTest((MonolithicKEpsilonHighReTest('testMonolithicKEpsilonHighReRfcVelocityTransient')))

    # adding monolithic k-omega tests
    nightSuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaAfcTkeSteady')))
    nightSuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaAfcVelocitySteady')))
    nightSuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaRfcTkeSteady')))
    nightSuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaRfcVelocitySteady')))
    nightSuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaRfcTkeTransient')))
    nightSuite.addTest((MonolithicKOmegaTest('testMonolithicKOmegaRfcVelocityTransient')))

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
