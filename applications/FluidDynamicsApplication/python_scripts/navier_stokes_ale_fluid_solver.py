from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication", "MeshMovingApplication")
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving

# other imports
from ale_fluid_solver import AleFluidSolver
import python_solvers_wrapper_fluid


def CreateSolver(model, solver_settings, parallelism):
    return NavierStokesAleFluidSolver(model, solver_settings, parallelism)


class NavierStokesAleFluidSolver(AleFluidSolver):
     def _CreateFluidSolver(self, solver_settings, parallelism):
        return python_solvers_wrapper_fluid.CreateSolverByParameters(
            self.model, solver_settings, parallelism)
