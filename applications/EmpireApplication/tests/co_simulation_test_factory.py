from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.EmpireApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from KratosMultiphysics.KratosUnittest import isclose as t_isclose
import os

import co_simulation_test_case

try:
    import numpy
    numpy_available = True
except ImportError:
    numpy_available = False

class TestSmallCoSimulationCases(co_simulation_test_case.CoSimulationTestCase):
    '''This class contains "small" CoSimulation-Cases, small enough to run in the nightly suite
    '''
    def test_MokFSI(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.createTest("fsi_mok", "cosim_mok_fsi")
            self.runTest()

class TestSmallCoSimulationPotentialCase(co_simulation_test_case.CoSimulationTestCase):
    '''This class contains "small" CoSimulation-Cases, small enough to run in the nightly suite
    '''
    def test_PotentialCase(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.createTest("fsi_potential_flow_sdof", "project_cosim_naca0012_small_fsi")
            self.runTestSteady()
            self._check_results(self.fluid_solver.lift_process.fluid_model_part.ProcessInfo.GetValue(CPFApp.LIFT_COEFFICIENT), 0.32782398801832946, 0.0, 1e-9)

class TestCoSimulationCases(co_simulation_test_case.CoSimulationTestCase):
    '''This class contains "full" CoSimulation-Cases, too large for the nightly suite and therefore
    have to be in the validation-suite
    '''
    def test_WallFSI(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.createTest("fsi_wall", "cosim_wall_weak_coupling_fsi")
            self.runTest()

    def test_SDoFDragRectangleFSI(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.createTest("fsi_sdof_drag_rectangle", "cosim_sdof_drag_rectangle_fsi")
            self.runTest()

    # def test_MDoFDragPitchRectangleFSI(self):
    #     with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
    #         self.createTest("fsi_mdof_drag_pitch_rectangle", "cosim_mdof_drag_pitch_rectangle_fsi")
    #         self.runTest()

    def _check_results(self, result, reference, rel_tol, abs_tol):
        isclosethis = t_isclose(result, reference, rel_tol, abs_tol)

        full_msg =  "Failed with following parameters:\n"
        full_msg += str(result) + " != " + str(reference) + ", rel_tol = "
        full_msg += str(rel_tol) + ", abs_tol = " + str(abs_tol)

        self.assertTrue(isclosethis, msg=full_msg)

if __name__ == '__main__':
    KratosUnittest.main()
