# import Kratos
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsHydraulicsApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from generalTests import KratosFluidDynamicsHydraulicsGeneralTests
from test_hydraulic_fluid_auxiliary_utilities import HydraulicFluidAuxiliaryUtilitiesTest
from apply_hydraulic_inlet_process_test import ApplyHydraulicInletProcessTest
from two_fluid_hydraulic_test import TwoFluidHydraulicSolverTest

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    # smallSuite will contain the following tests:
    # - testSmallExample
    smallSuite = suites['small']
    smallSuite.addTest(KratosFluidDynamicsHydraulicsGeneralTests('testSmallExample'))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ApplyHydraulicInletProcessTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([HydraulicFluidAuxiliaryUtilitiesTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TwoFluidHydraulicSolverTest]))


    # Create a test suit with the selected tests
    # nightSuite will contain the following tests:
    # - testSmallExample
    # - testNightlyFirstExample
    # - testNightlySecondExample
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            KratosFluidDynamicsHydraulicsGeneralTests
        ])
    )

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
