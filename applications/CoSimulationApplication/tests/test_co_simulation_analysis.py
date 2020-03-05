from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis
from unittest.mock import MagicMock


class TestCoSimulationAnalysis(KratosUnittest.TestCase):

    def test_GetSolver_one_level(self):
        co_sim_analysis = CoSimulationAnalysis(params)

        top_level_solver_wrapper = co_sim_analysis._GetSolver()

    def test_GetSolver_two_levels(self):
        co_sim_analysis = CoSimulationAnalysis(params)

        top_level_solver_wrapper  = co_sim_analysis._GetSolver()
        fluid_solver_wrapper      = co_sim_analysis._GetSolver("fluid")
        structural_solver_wrapper = co_sim_analysis._GetSolver("structure")

    def test_GetSolver_three_levels(self):
        co_sim_analysis = CoSimulationAnalysis(params)

        top_level_solver_wrapper  = co_sim_analysis._GetSolver()
        controller_solver_wrapper = co_sim_analysis._GetSolver("controller")
        fsi_solver_wrapper        = co_sim_analysis._GetSolver("fsi")
        fluid_solver_wrapper      = co_sim_analysis._GetSolver("fsi.fluid")
        structural_solver_wrapper = co_sim_analysis._GetSolver("fsi.structure")

        fsi_fluid_solver_wrapper      = fsi_solver_wrapper._GetSolver("fluid")
        fsi_structural_solver_wrapper = fsi_solver_wrapper._GetSolver("structure")


if __name__ == '__main__':
    KratosUnittest.main()
