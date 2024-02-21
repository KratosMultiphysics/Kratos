# import Kratos
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication         as MeshingApplication
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from test_refine import TestRedistance                                     as TTestRedistance
from test_remesh_rectangle import TestRemeshMMG2D                          as TTestRemeshMMG2D
from test_remesh_sphere import TestRemeshMMG3D                             as TTestRemeshMMG3D
from meshing_application_test_factory  import TwoDDynamicBeamTest          as TTwoDDynamicBeamTest
from meshing_application_test_factory  import TwoDDynamicBeamLineLoadTest  as TTwoDDynamicBeamLineLoadTest
from meshing_application_test_factory  import ThreeDShellTest              as TThreeDShellTest
from meshing_application_test_factory  import ThreeDDynamicBeamTest        as TThreeDDynamicBeamTest
from test_local_refine_parallel_to_boundaries import TestLocalRefineParallelToBoundaries as TTestRefineOnBoundaries
from test_local_refine_triangle_conditions import TestLocalRefineTriangleMeshConditions as TTestLocalRefineTriangleMeshConditions
from test_local_refine_only_on_boundaries import TestLocalRefineOnlyOnBoundaries as TTestLocalRefineOnlyOnBoundaries
from test_gradual_variable_interpolation_process import TestGradualVariableInterpolationProcess as TTestGradualVariableInterpolationProcess
from test_convert_linear_tetrahedra_to_quadratic_modeler import TestConvertLinearTetrahedraToQuadraticModeler
## NIGHTLY TESTS

## VALIDATION TESTS

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
    smallSuite = suites['small']
    smallSuite.addTest(TTestRefineOnBoundaries('test_refine_boundary_elems'))
    smallSuite.addTest(TTestLocalRefineTriangleMeshConditions('test_refine_condition_mesh'))
    smallSuite.addTest(TTestLocalRefineOnlyOnBoundaries('test_refine_on_boundary_edges'))
    smallSuite.addTest(TTestGradualVariableInterpolationProcess('test_gradual_variable_interpolation_process'))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([
        TestConvertLinearTetrahedraToQuadraticModeler]))
    if  hasattr(MeshingApplication,  "TetrahedraReconnectUtility") :
        smallSuite.addTest(TTestRedistance('test_refine_all'))
        smallSuite.addTest(TTestRedistance('test_refine_half'))
        smallSuite.addTest(TTestRedistance('test_refine_half_and_improve'))
    else:
        KratosMultiphysics.Logger.PrintWarning("Unittests", "TetrahedraReconnectUtility process is not compiled and the corresponding tests will not be executed")
    if hasattr(MeshingApplication,  "MmgProcess2D"):
        smallSuite.addTest(TTestRemeshMMG2D('test_remesh_rectangle_hessian'))
        smallSuite.addTest(TTestRemeshMMG3D('test_remesh_sphere'))
        smallSuite.addTest(TTestRemeshMMG3D('test_remesh_sphere_skin'))
        smallSuite.addTest(TTestRemeshMMG3D('test_remesh_sphere_skin_prisms'))
        smallSuite.addTest(TTestRemeshMMG3D('test_isosurface_remesh_sphere'))
        smallSuite.addTest(TTwoDDynamicBeamTest('test_execution'))
        smallSuite.addTest(TTwoDDynamicBeamLineLoadTest('test_execution'))
        smallSuite.addTest(TThreeDShellTest('test_execution'))
        smallSuite.addTest(TThreeDDynamicBeamTest('test_execution'))
    else:
        KratosMultiphysics.Logger.PrintWarning("Unittests", "MMG process is not compiled and the corresponding tests will not be executed")

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTest(TTestRefineOnBoundaries('test_refine_boundary_elems'))
    allSuite.addTest(TTestLocalRefineTriangleMeshConditions('test_refine_condition_mesh'))
    allSuite.addTest(TTestLocalRefineOnlyOnBoundaries('test_refine_on_boundary_edges'))
    allSuite.addTest(TTestGradualVariableInterpolationProcess('test_gradual_variable_interpolation_process'))
    if  hasattr(MeshingApplication, "TetrahedraReconnectUtility"):
        allSuite.addTests(
            KratosUnittest.TestLoader().loadTestsFromTestCases([
                TTestRedistance
            ])
        )
    else:
        KratosMultiphysics.Logger.PrintWarning("Unittests", "TetrahedraReconnectUtility process is not compiled and the corresponding tests will not be executed")

    if hasattr(MeshingApplication, "MmgProcess2D") :
        allSuite.addTests(
            KratosUnittest.TestLoader().loadTestsFromTestCases([
                TTestRemeshMMG2D,
                TTestRemeshMMG3D,
                TTwoDDynamicBeamTest,
                TTwoDDynamicBeamLineLoadTest,
                TThreeDShellTest,
                TThreeDDynamicBeamTest,
            ])
        )
    else:
        KratosMultiphysics.Logger.PrintWarning("Unittests", "MMG process is not compiled and the corresponding tests will not be executed")

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
