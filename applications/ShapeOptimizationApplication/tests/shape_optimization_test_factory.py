from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.kratos_utilities as kratos_utilities

# Other imports
import os

# Check if external Apps are available
has_eigen_app =  kratos_utilities.CheckIfApplicationsAvailable("EigenSolversApplication")
has_csm_app = kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication")
has_mesh_moving_app = kratos_utilities.CheckIfApplicationsAvailable("MeshMovingApplication")
has_mapping_app = kratos_utilities.CheckIfApplicationsAvailable("MappingApplication")

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

class mapper_test(ShapeOptimizationTestFactory):
    execution_directory = "mapper_test"
    execution_file = "run_test"

@kratos_unittest.skipUnless(has_csm_app,"Missing required application: StructuralMechanicsApplication")
class opt_process_shell_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_shell_test"
    execution_file = "run_test"

@kratos_unittest.skipUnless(has_csm_app and has_mesh_moving_app,"Missing (one or all) required applications: StructuralMechanicsApplication, MeshMovingApplication")
class opt_process_solid_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_solid_test"
    execution_file = "run_test"

@kratos_unittest.skipUnless(has_csm_app and has_eigen_app,"Missing (one or all) required applications: StructuralMechanicsApplication, EigenSolversApplication")
class opt_process_eigenfrequency_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_eigenfrequency_test"
    execution_file = "run_test"

@kratos_unittest.skipUnless(has_csm_app and has_eigen_app,"Missing (one or all) required applications: StructuralMechanicsApplication, EigenSolversApplication")
class opt_process_weighted_eigenfrequency_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_weighted_eigenfrequency_test"
    execution_file = "run_test"

@kratos_unittest.skipUnless(has_csm_app,"Missing required application: StructuralMechanicsApplication")
class algorithm_steepest_descent_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_steepest_descent_test"
    execution_file = "run_test"

@kratos_unittest.skipUnless(has_csm_app,"Missing required application: StructuralMechanicsApplication")
class algorithm_penalized_projection_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_penalized_projection_test"
    execution_file = "run_test"

@kratos_unittest.skipUnless(has_csm_app,"Missing required application: StructuralMechanicsApplication")
class algorithm_trust_region_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_trust_region_test"
    execution_file = "run_test"

class trust_region_projector_test(ShapeOptimizationTestFactory):
    execution_directory = "trust_region_projector_test"
    execution_file = "run_test"

@kratos_unittest.skipUnless(has_csm_app,"Missing required application: StructuralMechanicsApplication")
class algorithm_bead_optimization_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_bead_optimization_test"
    execution_file = "run_test"

@kratos_unittest.skipUnless(has_csm_app,"Missing required application: StructuralMechanicsApplication")
class opt_process_step_adaption_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_step_adaption_test"
    execution_file = "run_test"

@kratos_unittest.skipUnless(has_csm_app,"Missing required application: StructuralMechanicsApplication")
class opt_process_multiobjective_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_multiobjective_test"
    execution_file = "run_test"

@kratos_unittest.skipUnless(has_csm_app,"Missing required application: StructuralMechanicsApplication")
class opt_process_stress_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_stress_test"
    execution_file = "run_test"

@kratos_unittest.skipUnless(has_csm_app,"Missing required application: StructuralMechanicsApplication")
class sensitivity_verification_semi_analytic_process_test(ShapeOptimizationTestFactory):
    execution_directory = "sensitivity_verification_process_test"
    execution_file = "run_semi_analytic_step_size_verification"

@kratos_unittest.skipUnless(has_csm_app,"Missing required application: StructuralMechanicsApplication")
class sensitivity_verification_in_design_space_process_test(ShapeOptimizationTestFactory):
    execution_directory = "sensitivity_verification_process_test"
    execution_file = "run_sensitivity_verification_in_design_space"

@kratos_unittest.skipUnless(has_csm_app,"Missing required application: StructuralMechanicsApplication")
class sensitivity_verification_in_geometry_space_process_test(ShapeOptimizationTestFactory):
    execution_directory = "sensitivity_verification_process_test"
    execution_file = "run_sensitivity_verification_in_geometry_space"

@kratos_unittest.skipUnless(has_mapping_app,"Missing required application: MappingApplication")
class in_plane_opt_test(ShapeOptimizationTestFactory):
    execution_directory = "in_plane_opt_test"
    execution_file = "run_test"

class packaging_mesh_based_test(ShapeOptimizationTestFactory):
    execution_directory = "packaging_mesh_based_test"
    execution_file = "run_test"

class packaging_plane_based_test(ShapeOptimizationTestFactory):
    execution_directory = "packaging_plane_based_test"
    execution_file = "run_test"


# ==============================================================================