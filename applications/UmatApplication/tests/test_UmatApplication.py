# import Kratos
import KratosMultiphysics

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits:


from test_J2_plasticity import TestJ2Plasticity

def AssembleTestSuites():

    # Suites to run
    suites = KratosUnittest.KratosSuites

    # SMALL TESTS
    #small_suite = suites['small']

    # NIGHT TESTS
    night_suite = suites['nightly']
    night_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestJ2Plasticity]))


    # ALL TESTS
    all_suite = suites['all']

    all_suite.addTests(night_suite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
