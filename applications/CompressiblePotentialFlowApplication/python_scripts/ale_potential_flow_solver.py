from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Other imports
from  KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
have_mesh_moving = CheckIfApplicationsAvailable("MeshMovingApplication")
if have_mesh_moving:
    from KratosMultiphysics.MeshMovingApplication.ale_fluid_solver import AleFluidSolver
else:
    raise Exception("In importing the AlePotentialFlowSolver: The solver requires the MeshMovingApplication, but this application is not available.")

import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_solver as potential_flow_solver

def CreateSolver(model, solver_settings, parallelism):
    return AlePotentialFlowSolver(model, solver_settings, parallelism)


class AlePotentialFlowSolver(AleFluidSolver):
    def __init__(self, model, solver_settings, parallelism):
        super(AlePotentialFlowSolver, self).__init__(model, solver_settings, parallelism)
        self.fluid_solver.min_buffer_size = 1

    def _CreateFluidSolver(self, solver_settings, parallelism):
        return potential_flow_solver.CreateSolver(self.model, solver_settings)

    def SolveSolutionStep(self):
        is_converged = True

        for mesh_solver in self.mesh_motion_solvers:
            is_converged &= mesh_solver.SolveSolutionStep()

        is_converged &= self.fluid_solver.SolveSolutionStep()

        return is_converged

    @classmethod
    def _ManipulateFluidSolverSettingsForReactionsComputation(cls, fluid_solver_settings):
        if fluid_solver_settings.Has("compute_reactions"):
            if fluid_solver_settings["compute_reactions"].GetBool():
                fluid_solver_settings["compute_reactions"].SetBool(False)
                warn_msg  = '"compute_reactions" is switched on for the fluid-solver, switching it off!'
                KM.Logger.PrintWarning("::[AlePotentialFlowSolver]::", warn_msg)
        else:
            fluid_solver_settings.AddEmptyValue("compute_reactions").SetBool(False)
            info_msg = 'Setting "compute_reactions" to false for the fluid-solver'
            KM.Logger.PrintInfo("::[AlePotentialFlowSolver]::", info_msg)

    @classmethod
    def _ManipulateMeshMotionSolverSettingsForMeshVelocityComputation(cls, fluid_solver_settings, mesh_motion_solver_settings):
        if mesh_motion_solver_settings.Has("calculate_mesh_velocity"):
            if mesh_motion_solver_settings["calculate_mesh_velocity"].GetBool():
                mesh_motion_solver_settings["calculate_mesh_velocity"].SetBool(False)
                warn_msg  = '"calculate_mesh_velocity" is switched on for the mesh-motion-solver, switching it off!'
                KM.Logger.PrintWarning("::[AlePotentialFlowSolver]::", warn_msg)
        else:
            mesh_motion_solver_settings.AddEmptyValue("calculate_mesh_velocity").SetBool(False)
            info_msg = 'Setting "calculate_mesh_velocity" to false for the mesh-motion-solver'
            KM.Logger.PrintInfo("::[AlePotentialFlowSolver]::", info_msg)