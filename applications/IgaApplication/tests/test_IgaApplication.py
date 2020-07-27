# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import Iga test factory tests
from iga_test_factory import SinglePatchTest as SinglePatchTest
from iga_test_factory import Shell3pLinearBeamThick as TShell3pLinearBeamThick
from iga_test_factory import Shell3pNonLinearBeamThick as TShell3pNonLinearBeamThick
from iga_test_factory import Shell3pNonLinearBeamThickSD as TShell3pNonLinearBeamThickSD
from iga_test_factory import Shell3pLinearScordelis as TShell3pLinearScordelis

# Import the tests o test_classes to create the suits

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

    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([
        SinglePatchTest,
        TShell3pLinearBeamThick,
        TShell3pNonLinearBeamThick,
        TShell3pNonLinearBeamThickSD,
        TShell3pLinearScordelis
        ]))


    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
