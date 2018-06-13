# import Kratos
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import KratosMultiphysics.PfemApplication as PfemApplication
import KratosMultiphysics.PfemFluidDynamicsApplication as PfemFluidDynamicsApplication     
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication

# Import the tests o test_classes to create the suits:

# NIGTHLY TESTS
import NightTests

def AssambleTestSuites():

    # Suites to run
    suites = KratosUnittest.KratosSuites

    # NIGHTLY TESTS
    night_suite = NightTests.SetTestSuite(suites)

    # ALL TESTS
    all_suite = suites['all']

    all_suite.addTests(night_suite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())

