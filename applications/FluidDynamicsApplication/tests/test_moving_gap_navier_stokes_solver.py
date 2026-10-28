#     |  /           |
#     ' /   __| _` | __|  _ \   __|
#     . \  |   (   | |   (   |\__ `
#    _|\_\_|  \__,_|\__|\___/ ____/
#                    Multi-Physics
#
#   Fluid Dynamics Application
#
#   License:         BSD License
#                    Kratos default license: kratos/license.txt
#
#   Main author:     Nicolò Antonelli
#

"""Configuration checks for the QSVMS Picard strategy."""
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication.moving_gap_navier_stokes_solver import (
    MovingGapNavierStokesSolver, create_qsvms_picard_strategy)


class MovingGapNavierStokesSolverTest(KratosUnittest.TestCase):
    def test_picard_reassembles_and_rejects_exact_history(self):
        model = KM.Model()
        part = model.CreateModelPart('Fluid')
        solver = MovingGapNavierStokesSolver(model, KM.Parameters('{}'))
        with self.assertRaisesRegex(ValueError, 'numerical history'):
            solver.SetExactHistoryCallback(lambda node, time: [0.0, 0.0, 0.0])
        self.assertEqual(solver.settings['coupling_settings']['formulation'].GetString(), 'qsvms')
        strategy = create_qsvms_picard_strategy(part, KM.SkylineLUFactorizationSolver(), solver.settings)
        self.assertEqual(strategy.GetRebuildLevel(), 2)
        strategy.Clear()


if __name__ == '__main__':
    KratosUnittest.main()
