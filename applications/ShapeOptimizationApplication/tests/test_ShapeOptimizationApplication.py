# ==============================================================================
# Imports
# ==============================================================================


# Import Kratos "wrapper" for unittests
import KratosMultiphysics as km
import KratosMultiphysics.KratosUnittest as KratosUnittest

# ==============================================================================
# Import the tests or test_classes to create the suits
# ==============================================================================

# Small tests
from shape_optimization_test_factory import opt_process_shell_test
from shape_optimization_test_factory import opt_process_solid_test
from shape_optimization_test_factory import opt_process_vertex_morphing_test
from shape_optimization_test_factory import opt_process_vertex_morphing_small_test
from shape_optimization_test_factory import opt_process_eigenfrequency_test
from shape_optimization_test_factory import opt_process_weighted_eigenfrequency_test
from shape_optimization_test_factory import algorithm_steepest_descent_test
from shape_optimization_test_factory import algorithm_penalized_projection_test
from shape_optimization_test_factory import algorithm_trust_region_test
from shape_optimization_test_factory import trust_region_projector_test
from shape_optimization_test_factory import algorithm_gradient_projection_test
from shape_optimization_test_factory import algorithm_qn_bb_relaxed_gradient_projection_test
from shape_optimization_test_factory import algorithm_bead_optimization_test
from shape_optimization_test_factory import algorithm_shape_fraction_test
from shape_optimization_test_factory import opt_process_step_adaption_test
from shape_optimization_test_factory import mapper_test
from shape_optimization_test_factory import opt_process_multiobjective_test
# from shape_optimization_test_factory import opt_process_stress_test
# from shape_optimization_test_factory import sensitivity_verification_semi_analytic_process_test
# from shape_optimization_test_factory import sensitivity_verification_in_design_space_process_test
# from shape_optimization_test_factory import sensitivity_verification_in_geometry_space_process_test
from shape_optimization_test_factory import in_plane_opt_test
from shape_optimization_test_factory import packaging_mesh_based_test
from shape_optimization_test_factory import packaging_plane_based_test
from shape_optimization_test_factory import remeshing_opt_process_test
from shape_optimization_test_factory import sliding_opt_test
from shape_optimization_test_factory import curvature_3NTriangle_test, curvature_6NTriangle_test, curvature_4NQuad_test, curvature_8NQuad_test
from shape_optimization_test_factory import mapper_adaptive_filter_curvature_test
from shape_optimization_test_factory import sensitivity_heatmap_test
from wrl_io_test.test_wrl_io import WrlIOTest
from surface_normal_shape_change_response_test.test_surface_normal_shape_change_response import SurfaceNormalShapeChangeTest
from face_angle_response_test.test_face_angle_response import FaceAngleTest
from mapper_plane_symmetry_test.plane_symmetry_test import PlaneSymmetryMapperTest
from mapper_revolution_test.revolution_test import RevolutionMapperTest
from shape_optimization_test_factory import direction_damping_test
from total_volume.test_total_volume_response import TestTotalVolumeResponseFunction2D
from total_volume.test_total_volume_response import TestTotalVolumeResponseFunction3D

# Nightly tests

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
    smallSuite.addTest(mapper_test('test_execution'))
    smallSuite.addTest(opt_process_vertex_morphing_small_test('test_execution'))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([WrlIOTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([SurfaceNormalShapeChangeTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FaceAngleTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([PlaneSymmetryMapperTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([RevolutionMapperTest]))
    smallSuite.addTest(algorithm_gradient_projection_test('test_execution'))
    smallSuite.addTest(algorithm_qn_bb_relaxed_gradient_projection_test('test_execution'))
    smallSuite.addTest(algorithm_shape_fraction_test('test_execution'))
    smallSuite.addTest(remeshing_opt_process_test('test_execution'))
    smallSuite.addTest(sliding_opt_test('test_execution'))
    smallSuite.addTest(direction_damping_test('test_execution'))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTotalVolumeResponseFunction2D]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTotalVolumeResponseFunction3D]))
    smallSuite.addTest(curvature_3NTriangle_test('test_execution'))
    smallSuite.addTest(curvature_6NTriangle_test('test_execution'))
    smallSuite.addTest(curvature_4NQuad_test('test_execution'))
    smallSuite.addTest(curvature_8NQuad_test('test_execution'))
    smallSuite.addTest(mapper_adaptive_filter_curvature_test('test_execution'))
    smallSuite.addTest(sensitivity_heatmap_test('test_execution'))

    # Adding nightly tests (tests that take < 10min)
    nightSuite = suites['nightly']
    nightSuite.addTest(opt_process_shell_test('test_execution'))
    nightSuite.addTest(opt_process_solid_test('test_execution'))
    nightSuite.addTest(algorithm_bead_optimization_test('test_execution'))
    nightSuite.addTest(opt_process_step_adaption_test('test_execution'))
    nightSuite.addTest(in_plane_opt_test('test_execution'))
    nightSuite.addTest(packaging_mesh_based_test('test_execution'))
    nightSuite.addTest(packaging_plane_based_test('test_execution'))
    nightSuite.addTest(opt_process_vertex_morphing_test('test_execution'))
    nightSuite.addTest(opt_process_eigenfrequency_test('test_execution'))
    nightSuite.addTest(opt_process_weighted_eigenfrequency_test('test_execution'))
    nightSuite.addTest(trust_region_projector_test('test_execution'))
    nightSuite.addTest(opt_process_multiobjective_test('test_execution'))
    # nightSuite.addTest(opt_process_stress_test('test_execution'))
    # nightSuite.addTest(sensitivity_verification_semi_analytic_process_test('test_execution'))
    # nightSuite.addTest(sensitivity_verification_in_design_space_process_test('test_execution'))
    # nightSuite.addTest(sensitivity_verification_in_geometry_space_process_test('test_execution'))

    # Adding small tests to nightly tests
    nightSuite.addTests(smallSuite)

    # Adding validation tests
    validationSuite = suites['validation']
    validationSuite.addTests(nightSuite)
    validationSuite.addTest(algorithm_trust_region_test('test_execution'))
    validationSuite.addTest(algorithm_steepest_descent_test('test_execution'))
    validationSuite.addTest(algorithm_penalized_projection_test('test_execution'))

    # Creating a test suit that contains all tests:
    allSuite = suites['all']
    allSuite.addTests(validationSuite)

    return suites

# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())

# ==============================================================================
