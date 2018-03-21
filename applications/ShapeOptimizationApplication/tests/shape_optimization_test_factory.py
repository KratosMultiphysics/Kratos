from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

# Other imports
import os
import sys
import shutil

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
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))+"/"+self.execution_directory):
            test_status = os.system("runkratos "+self.execution_file+" >run_test.log")
            if test_status == 0:
                kratos_utils.DeleteFileIfExisting("run_test.log")
            else:
                raise RuntimeError("test_run.py failed!")

    # --------------------------------------------------------------------------
    def tearDown(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            kratos_utils.DeleteDirectoryIfExisting("__pycache__")
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))+"/"+self.execution_directory):
            kratos_utils.DeleteDirectoryIfExisting("__pycache__")

# ==============================================================================
class opt_process_vertex_morphing_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_vertex_morphing_test"
    execution_file = "run_test.py"

class opt_process_shell_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_shell_test"
    execution_file = "run_test.py"

class opt_process_solid_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_solid_test"
    execution_file = "run_test.py"

class opt_process_eigenfrequency_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_eigenfrequency_test"
    execution_file = "run_test.py"

class algorithm_steepest_descent_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_steepest_descent_test"
    execution_file = "run_test.py"

class algorithm_penalized_projection_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_penalized_projection_test"
    execution_file = "run_test.py"

# ==============================================================================