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

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # adding incompressible potential flow solver tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([IncompressiblePotentialFlowSolverFormulationTest]))

    # adding monolithic flow solver tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicVelocityPressureFormulationTest]))

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
