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
from KratosMultiphysics.FluidDynamicsApplication.moving_gap_lifecycle import MovingGapLifecycle


class MovingGapNavierStokesSolverTest(KratosUnittest.TestCase):
    def test_lifecycle_rejects_invalid_step_before_mutating_model(self):
        lifecycle = MovingGapLifecycle.__new__(MovingGapLifecycle)
        lifecycle._pending = False
        lifecycle._failed = False
        lifecycle._last_time = 3.0
        lifecycle.dt = 0.1
        with self.assertRaisesRegex(ValueError, 'fixed time steps'):
            lifecycle.prepare_solution_step(3.2)
        self.assertFalse(lifecycle._failed)
        lifecycle._failed = True
        with self.assertRaisesRegex(RuntimeError, 'restore a checkpoint'):
            lifecycle.prepare_solution_step(3.1)

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
