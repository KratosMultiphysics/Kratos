from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os
import co_simulation_test_case


class TestCoSimulationNightlyCases(co_simulation_test_case.CoSimulationTestCase):
    # This class contains "nightly" CoSimulation cases, small enough to run in the nightly build.
    def test_MokFSI(self):
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.createTest("fsi_mok", "cosim_mok_fsi")
            self.runTest()


class TestCoSimulationValidationCases(co_simulation_test_case.CoSimulationTestCase):
    # This class contains "validation" CoSimulation cases, too large for the nightly suite and therefore have to be in
    # the validation suite.
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
