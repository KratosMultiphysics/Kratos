from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Other imports
from KratosMultiphysics.MeshMovingApplication.ale_fluid_solver import AleFluidSolver
import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_solver as potential_flow_solver

def CreateSolver(model, solver_settings, parallelism):
    return AlePotentialFlowSolver(model, solver_settings, parallelism)


class AlePotentialFlowSolver(AleFluidSolver):
    def __init__(self, model, solver_settings, parallelism):
        super(AlePotentialFlowSolver, self).__init__(model, solver_settings, parallelism)
        self.fluid_solver.min_buffer_size = 1
        self.settings["fluid_solver_settings"]["compute_reactions"].SetBool(False)

    def _CreateFluidSolver(self, solver_settings, parallelism):
        return potential_flow_solver.CreateSolver(self.model, solver_settings)

    def SolveSolutionStep(self):
        is_converged = True

        for mesh_solver in self.mesh_motion_solvers:
            is_converged &= mesh_solver.SolveSolutionStep()

        is_converged &= self.fluid_solver.SolveSolutionStep()

        return is_converged
