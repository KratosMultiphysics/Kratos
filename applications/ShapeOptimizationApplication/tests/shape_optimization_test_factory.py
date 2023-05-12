
import KratosMultiphysics as KM

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.kratos_utilities as kratos_utilities

# Other imports
import os

try:
    import KratosMultiphysics.MeshingApplication
    has_mmg = hasattr(KratosMultiphysics.MeshingApplication, "MmgProcess2D")
except ImportError:
    has_mmg = False

# ==============================================================================
class ShapeOptimizationTestFactory(kratos_unittest.TestCase):
    # --------------------------------------------------------------------------
    def setUp(self):
        KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

    # --------------------------------------------------------------------------
    def test_execution(self):
        with kratos_unittest.WorkFolderScope(self.execution_directory, __file__):
            __import__(self.execution_directory+"."+self.execution_file)

    # --------------------------------------------------------------------------
    def tearDown(self):
        with kratos_unittest.WorkFolderScope(self.execution_directory, __file__):
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")

# ==============================================================================
class opt_process_vertex_morphing_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_vertex_morphing_test"
    execution_file = "run_test"

class opt_process_vertex_morphing_small_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_vertex_morphing_test"
    execution_file = "run_test_small"

class mapper_test(ShapeOptimizationTestFactory):
    execution_directory = "mapper_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class opt_process_shell_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_shell_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "MeshMovingApplication")
class opt_process_solid_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_solid_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "LinearSolversApplication")
class opt_process_eigenfrequency_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_eigenfrequency_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "LinearSolversApplication")
class opt_process_weighted_eigenfrequency_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_weighted_eigenfrequency_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class algorithm_steepest_descent_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_steepest_descent_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class algorithm_penalized_projection_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_penalized_projection_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class algorithm_gradient_projection_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_gradient_projection_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class algorithm_qn_bb_relaxed_gradient_projection_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_qn_bb_relaxed_gradient_projection_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class algorithm_trust_region_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_trust_region_test"
    execution_file = "run_test"

class trust_region_projector_test(ShapeOptimizationTestFactory):
    execution_directory = "trust_region_projector_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class algorithm_bead_optimization_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_bead_optimization_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class algorithm_shape_fraction_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_shape_fraction_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class opt_process_step_adaption_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_step_adaption_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class opt_process_multiobjective_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_multiobjective_test"
    execution_file = "run_test"

# @kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
# class opt_process_stress_test(ShapeOptimizationTestFactory):
#     execution_directory = "opt_process_stress_test"
#     execution_file = "run_test"

# @kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
# class sensitivity_verification_semi_analytic_process_test(ShapeOptimizationTestFactory):
#     execution_directory = "sensitivity_verification_process_test"
#     execution_file = "run_semi_analytic_step_size_verification"

# @kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
# class sensitivity_verification_in_design_space_process_test(ShapeOptimizationTestFactory):
#     execution_directory = "sensitivity_verification_process_test"
#     execution_file = "run_sensitivity_verification_in_design_space"

# @kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
# class sensitivity_verification_in_geometry_space_process_test(ShapeOptimizationTestFactory):
#     execution_directory = "sensitivity_verification_process_test"
#     execution_file = "run_sensitivity_verification_in_geometry_space"

@kratos_unittest.skipIfApplicationsNotAvailable("MappingApplication")
class in_plane_opt_test(ShapeOptimizationTestFactory):
    execution_directory = "in_plane_opt_test"
    execution_file = "run_test"

class packaging_mesh_based_test(ShapeOptimizationTestFactory):
    execution_directory = "packaging_mesh_based_test"
    execution_file = "run_test"

class packaging_plane_based_test(ShapeOptimizationTestFactory):
    execution_directory = "packaging_plane_based_test"
    execution_file = "run_test"

@kratos_unittest.skipUnless(has_mmg, "Test requires mmg library")
class remeshing_opt_process_test(ShapeOptimizationTestFactory):
    execution_directory = "remeshing_opt_process_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "MappingApplication")
class sliding_opt_test(ShapeOptimizationTestFactory):
    execution_directory = "sliding_opt_test"
    execution_file = "run_test"

class direction_damping_test(ShapeOptimizationTestFactory):
    execution_directory = "direction_damping_test"
    execution_file = "run_test"

class curvature_3NTriangle_test(ShapeOptimizationTestFactory):
    execution_directory = "geom_util_curvature_3NTriangle_test"
    execution_file = "run_test"

class curvature_6NTriangle_test(ShapeOptimizationTestFactory):
    execution_directory = "geom_util_curvature_6NTriangle_test"
    execution_file = "run_test"

class curvature_4NQuad_test(ShapeOptimizationTestFactory):
    execution_directory = "geom_util_curvature_4NQuad_test"
    execution_file = "run_test"

class curvature_8NQuad_test(ShapeOptimizationTestFactory):
    execution_directory = "geom_util_curvature_8NQuad_test"
    execution_file = "run_test"

class mapper_adaptive_filter_curvature_test(ShapeOptimizationTestFactory):
    execution_directory = "mapper_adaptive_filter_curvature_test"
    execution_file = "run_test"

class sensitivity_heatmap_test(ShapeOptimizationTestFactory):
    execution_directory = "sensitivity_heatmap_test"
    execution_file = "run_test"
# ==============================================================================