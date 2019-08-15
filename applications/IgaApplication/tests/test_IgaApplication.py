# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from node_curve_geometry_3d_tests import NodeCurveGeometry3DTests
from node_surface_geometry_3d_tests import NodeSurfaceGeometry3DTests
from iga_truss_element_tests import IgaTrussElementTests
from shell_kl_discrete_element_tests import ShellKLDiscreteElementTests
#from test_iga_output_process import TestIgaOutputProcess
# from test_iga_test_factory import IgaTestFactory
from iga_3p_element_tests import Iga3pElementTests
from iga_5p_element_tests import Iga5pElementTests

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
        NodeCurveGeometry3DTests,
        NodeSurfaceGeometry3DTests,
        IgaTrussElementTests,
        ShellKLDiscreteElementTests,
#        TestIgaOutputProcess,
        Iga3pElementTests,
        Iga5pElementTests
    ]))

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    #run_cpp_unit_tests.run()
    KratosUnittest.runTests(AssembleTestSuites())
