from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities

# Other imports
import os

class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

# ==============================================================================
class ShapeOptimizationTestFactory(KratosUnittest.TestCase):
    # --------------------------------------------------------------------------
    def setUp(self):
        pass

    # --------------------------------------------------------------------------
    def test_execution(self):
        test_scope = os.path.join(os.path.dirname(os.path.realpath(__file__)), self.execution_directory)
        with controlledExecutionScope(test_scope):
            __import__(self.execution_directory+"."+self.execution_file)

    # --------------------------------------------------------------------------
    def tearDown(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")

        test_scope = os.path.join(os.path.dirname(os.path.realpath(__file__)), self.execution_directory)
        with controlledExecutionScope(test_scope):
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")

# ==============================================================================
class opt_process_vertex_morphing_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_vertex_morphing_test"
    execution_file = "run_test"

class opt_process_shell_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_shell_test"
    execution_file = "run_test"

class opt_process_solid_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_solid_test"
    execution_file = "run_test"

class opt_process_eigenfrequency_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_eigenfrequency_test"
    execution_file = "run_test"

class opt_process_weighted_eigenfrequency_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_weighted_eigenfrequency_test"
    execution_file = "run_test"

class algorithm_steepest_descent_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_steepest_descent_test"
    execution_file = "run_test"

class algorithm_penalized_projection_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_penalized_projection_test"
    execution_file = "run_test"

class algorithm_trust_region_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_trust_region_test"
    execution_file = "run_test"

class trust_region_projector_test(ShapeOptimizationTestFactory):
    execution_directory = "trust_region_projector_test"
    execution_file = "run_test"

class mapper_test(ShapeOptimizationTestFactory):
    execution_directory = "mapper_test"
    execution_file = "run_test"

# ==============================================================================