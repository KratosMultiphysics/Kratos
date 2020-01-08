# import Kratos
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import KratosMultiphysics.DelaunayMeshingApplication as DelaunayMeshingApplication
import KratosMultiphysics.PfemFluidDynamicsApplication as PfemFluidDynamicsApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits:

# SMALL TESTS
import SmallTests

# NIGTHLY TESTS
import NightTests

def AssembleTestSuites():

    # Suites to run
    suites = KratosUnittest.KratosSuites

    # SMALL TESTS
    small_suite = SmallTests.SetTestSuite(suites)

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

