from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.EmpireApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from compare_two_files_check_process import CompareTwoFilesCheckProcess

import os

import co_simulation_test_case

try:
    import scipy
    import sympy
    scipy_and_sympy_available = True
except ImportError:
    scipy_and_sympy_available = False

try:
    import numpy
    numpy_available = True
except ImportError:
    numpy_available = False


def compareResults(reference_file, results_file):
    settings_check_process = KratosMultiphysics.Parameters("""
    {
        "reference_file_name"   : "",
        "output_file_name"      : "",
        "comparison_type"       : "dat_file",
        "remove_output_file"    : true,
        "tolerance"             : 1e-6
    }
    """)

    settings_check_process["reference_file_name"].SetString(reference_file)
    settings_check_process["output_file_name"].SetString(results_file)

    # creating a dummy model
    check_process = CompareTwoFilesCheckProcess(settings_check_process)

    check_process.ExecuteInitialize()
    check_process.ExecuteBeforeSolutionLoop()
    check_process.ExecuteInitializeSolutionStep()
    check_process.ExecuteFinalizeSolutionStep()
    check_process.ExecuteBeforeOutputStep()
    check_process.ExecuteAfterOutputStep()
    check_process.ExecuteFinalize()

class TestKratosSolver(co_simulation_test_case.CoSimulationTestCase):
    def test_KratosStructuralMechanicsSolver(self):
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # self.createTest('test_structural_mesh_motion_2d/rectangle_2D3N_test')
            # self.runTest()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D3N_test.time")

    def test_KratosFluidDynamicsSolver(self):
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # self.createTest('test_structural_mesh_motion_2d/rectangle_2D3N_test')
            # self.runTest()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D3N_test.time")

class TestSDoFSolver(co_simulation_test_case.CoSimulationTestCase):
    def test_SDoFSolver(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            folder_name = "sdof_solver"
            self.createTest("sdof_solver", "cosim_sdof")
            self.runTest()
            reference_file = os.path.join(folder_name,"results_sdof_ref.dat")
            result_file = os.path.join(folder_name,"results_sdof.dat")
            compareResults(reference_file, result_file)

class TestMDoFSolver(co_simulation_test_case.CoSimulationTestCase):
    def test_MDoFSDoFModel(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            folder_name = "mdof_solver"
            self.createTest(folder_name, "cosim_mdof_sdof")
            self.runTest()
            reference_file = os.path.join(folder_name,"results_mdof_sdof_ref.dat")
            result_file = os.path.join(folder_name,"results_mdof_sdof.dat")
            compareResults(reference_file, result_file)

    def test_MDoFGenericModel(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            folder_name = "mdof_solver"
            self.createTest(folder_name, "cosim_mdof_generic")
            self.runTest()
            reference_file = os.path.join(folder_name,"results_mdof_generic_ref.dat")
            result_file = os.path.join(folder_name,"results_mdof_generic.dat")
            compareResults(reference_file, result_file)

    def test_MDoFCantileverShear2DModel(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not scipy_and_sympy_available:
            self.skipTest("Scipy/Sympy not available")
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            folder_name = "mdof_solver"
            self.createTest(folder_name, "cosim_mdof_cantilever_shear_2d")
            self.runTest()
            reference_file = os.path.join(folder_name,"results_mdof_cantilever_shear_2d_ref.dat")
            result_file = os.path.join(folder_name,"results_mdof_cantilever_shear_2d.dat")
            compareResults(reference_file, result_file)

    def test_MDoFBridge2DoFModel(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not scipy_and_sympy_available:
            self.skipTest("Scipy/Sympy not available")
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            folder_name = "mdof_solver"
            self.createTest(folder_name, "cosim_mdof_bridge_2dof")
            self.runTest()
            reference_file = os.path.join(folder_name,"results_mdof_bridge_2dof_ref.dat")
            result_file = os.path.join(folder_name,"results_mdof_bridge_2dof.dat")
            compareResults(reference_file, result_file)

class TestEmpireSolver(co_simulation_test_case.CoSimulationTestCase):
    def test_EmpireSolverWrapper(self):
        if "EMPIRE_API_LIBSO_ON_MACHINE" not in os.environ:
            self.skipTest("EMPIRE is not available")
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # self.createTest('test_structural_mesh_motion_2d/rectangle_2D3N_test')
            # self.runTest()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D3N_test.time")


if __name__ == '__main__':
    KratosUnittest.main()
