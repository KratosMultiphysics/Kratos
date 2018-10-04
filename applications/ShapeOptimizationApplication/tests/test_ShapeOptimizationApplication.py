# ==============================================================================
# Imports
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics import ShapeOptimizationApplication
from KratosMultiphysics import StructuralMechanicsApplication
from KratosMultiphysics import ExternalSolversApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import necessary external applications
try:
    from KratosMultiphysics import EigenSolversApplication
    is_eigen_app_missing = False
except ImportError as e:
    print("WARNING: EigenSolversApplication is not available, skip related tests!")
    is_eigen_app_missing = True
    
try:
    from KratosMultiphysics import MeshMovingApplication
    is_mesh_moving_app_missing = False
except ImportError as e:
    print("WARNING: MeshMovingApplication is not available, skip related tests!")
    is_mesh_moving_app_missing = True

# ==============================================================================
# Import the tests or test_classes to create the suits
# ==============================================================================

# Small tests
from shape_optimization_test_factory import opt_process_shell_test as opt_process_shell_test
from shape_optimization_test_factory import opt_process_solid_test as opt_process_solid_test
from shape_optimization_test_factory import opt_process_vertex_morphing_test as opt_process_vertex_morphing_test
from shape_optimization_test_factory import opt_process_eigenfrequency_test as opt_process_eigenfrequency_test
from shape_optimization_test_factory import opt_process_weighted_eigenfrequency_test as opt_process_weighted_eigenfrequency_test
from shape_optimization_test_factory import algorithm_steepest_descent_test as algorithm_steepest_descent_test
from shape_optimization_test_factory import algorithm_penalized_projection_test as algorithm_penalized_projection_test
from shape_optimization_test_factory import algorithm_trust_region_test as algorithm_trust_region_test
from shape_optimization_test_factory import trust_region_projector_test as trust_region_projector_test

# Niglty tests

# Validation tests

# ==============================================================================
# Test assembly
# ==============================================================================
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

    # Adding small tests (tests that take < 1s)
    smallSuite = suites['small']
    smallSuite.addTest(opt_process_vertex_morphing_test('test_execution'))
    smallSuite.addTest(opt_process_shell_test('test_execution'))
    if not is_mesh_moving_app_missing:
        smallSuite.addTest(opt_process_solid_test('test_execution'))
    if not is_eigen_app_missing:
        smallSuite.addTest(opt_process_eigenfrequency_test('test_execution'))
        smallSuite.addTest(opt_process_weighted_eigenfrequency_test('test_execution'))
    smallSuite.addTest(algorithm_steepest_descent_test('test_execution'))
    smallSuite.addTest(algorithm_penalized_projection_test('test_execution'))
    smallSuite.addTest(algorithm_trust_region_test('test_execution'))
    smallSuite.addTest(trust_region_projector_test('test_execution'))

    # Adding nightly tests (tests that take < 10min)
    nightSuite = suites['nightly']

    # Adding small tests to nightly tests
    nightSuite.addTests(smallSuite)

    # Adding validation tests
    validationSuite = suites['validation']

    # Creating a test suit that contains all tests:
    allSuite = suites['all']
    # allSuite.addTests(smallSuite) #Already added to small tests
    allSuite.addTests(nightSuite)
    allSuite.addTests(validationSuite)

    return suites

# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())

# ==============================================================================