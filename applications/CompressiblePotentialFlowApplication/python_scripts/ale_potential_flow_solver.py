from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# other imports
from KratosMultiphysics.MeshMovingApplication.ale_fluid_solver import AleFluidSolver
import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_solver as potential_flow_solver

def CreateSolver(model, solver_settings, parallelism):
    return AlePotentialFlowSolver(model, solver_settings, parallelism)


class AlePotentialFlowSolver(AleFluidSolver):
    def _CreateFluidSolver(self, solver_settings, parallelism):
        return potential_flow_solver.CreateSolver(self.model, solver_settings)