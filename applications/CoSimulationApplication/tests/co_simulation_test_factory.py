from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.CoSimulationApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

import os

import co_simulation_test_case

class TestSmallCoSimulationCases(co_simulation_test_case.CoSimulationTestCase):
    '''This class contains "small" CoSimulation-Cases, small enough to run in the nightly suite
    '''
    def test_MokFSI(self):
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.createTest("fsi_mok", "cosim_mok_fsi")
            self.runTest()

class TestCoSimulationCases(co_simulation_test_case.CoSimulationTestCase):
    '''This class contains "full" CoSimulation-Cases, too large for the nightly suite and therefore
    have to be in the validation-suite
    '''
    def test_WallFSI(self):
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.createTest("fsi_wall", "cosim_wall_weak_coupling_fsi")
            self.runTest()

    def test_SDoFDragRectangleFSI(self):
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.createTest("fsi_sdof_drag_rectangle", "cosim_sdof_drag_rectangle_fsi")
            self.runTest()

if __name__ == '__main__':
    KratosUnittest.main()
