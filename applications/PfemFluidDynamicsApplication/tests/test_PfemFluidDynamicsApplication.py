# import Kratos
import KratosMultiphysics

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suites:

# SMALL TESTS
import SmallTests

# NIGTHLY TESTS
import NightTests
from test_cut_pfem_2d import TestCutPfem

def AssembleTestSuites():

    # Suites to run
    suites = KratosUnittest.KratosSuites

    # SMALL TESTS
    small_suite = SmallTests.SetTestSuite(suites)
    small_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCutPfem]))

    # NIGHTLY TESTS
    night_suite = NightTests.SetTestSuite(suites)

    # inlude small suite in night suite
    night_suite.addTests(small_suite)

    # ALL TESTS
    all_suite = suites['all']

    all_suite.addTests(night_suite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())

