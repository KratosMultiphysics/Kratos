# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import Iga test factory tests
from iga_test_factory import SinglePatchTest as SinglePatchTest

# Import python based Iga simulation test
from shell_3p_element_tests import Shell3pElementTests

# Import the tests o test_classes to create the suits

def AssembleTestSuites():
    '''
    This application provides some tests for the element formulation
    of this application. It is recommended to use the IgaTestfactory
    to add more tests.
    It can be necessary to have the StructuralMechanicsApplication
    compiled to run all tests.
    The IgaApplication is part of the CI, thus tests in the "small"-suite
    need to run very fast.
    '''

    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([
        SinglePatchTest,
        Shell3pElementTests
        ]))


    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
