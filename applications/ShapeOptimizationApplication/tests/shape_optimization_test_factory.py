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

# ==============================================================================
class ShapeOptimizationTestFactory(KratosUnittest.TestCase):
    # --------------------------------------------------------------------------
    def setUp(self):
        pass

    # --------------------------------------------------------------------------
    def test_execution(self):
        original_directory = os.getcwd()
        os.chdir(self.execution_directory)

        __import__(self.execution_directory+"."+self.execution_file)

        os.chdir(original_directory)

    # --------------------------------------------------------------------------
    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting("__pycache__")

# ==============================================================================
class algorithm_steepest_descent_test(ShapeOptimizationTestFactory):
    execution_directory = "algorithm_steepest_descent_test"
    execution_file = "run_test"

class opt_process_shell_test(ShapeOptimizationTestFactory):
    execution_directory = "opt_process_shell_test"
    execution_file = "run_test"

# ==============================================================================