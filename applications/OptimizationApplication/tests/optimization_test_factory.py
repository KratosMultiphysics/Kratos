
import KratosMultiphysics as KM

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.kratos_utilities as kratos_utilities

# Other imports

# ==============================================================================
class OptimizationTestFactory(kratos_unittest.TestCase):
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

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "LinearSolversApplication")
class top_opt_test(OptimizationTestFactory):
    execution_directory = "top_opt_test"
    execution_file = "run_test"

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "LinearSolversApplication")
class mat_opt_test(OptimizationTestFactory):
    execution_directory = "mat_opt_test"
    execution_file = "run_test"    
    
@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "LinearSolversApplication")
class shell_shape_opt_test(OptimizationTestFactory):
    execution_directory = "shell-shape-opt-test"
    execution_file = "run_test"   

@kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "LinearSolversApplication")
class shell_thick_opt_test(OptimizationTestFactory):
    execution_directory = "shell-thickness-opt-test"
    execution_file = "run_test"        

# ==============================================================================
