from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MeshMovingApplication", "FluidDynamicsApplication","ExternalSolversApplication")

# Import applications
import KratosMultiphysics.MeshMovingApplication as MeshMovingApplication
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication

# Other imports
import navier_stokes_solver_fractionalstep
#import mesh_solver_base


def CreateSolver(model_part, custom_settings):
    return ALENavierStokesSolverFractionalStep(model_part, custom_settings)


class ALENavierStokesSolverFractionalStep(navier_stokes_solver_fractionalstep.NavierStokesSolverFractionalStep):
    def __init__(self, model_part, custom_settings):
        # remove the ale_settings so we can use the navier stokes constructor
        navier_stokes_settings = custom_settings.Clone()
        navier_stokes_settings.RemoveValue("ale_settings")
        navier_stokes_settings["solver_type"].SetString("FractionalStep")
        super(ALENavierStokesSolverFractionalStep, self).__init__(model_part, navier_stokes_settings)
        # create mesh motion solver
        custom_settings.AddEmptyValue("problem_data")
        custom_settings["problem_data"].AddEmptyValue("parallel_type")
        custom_settings["problem_data"]["parallel_type"].SetString("OpenMP")
        custom_settings.AddValue("solver_settings", custom_settings["ale_settings"])
        custom_settings.RemoveValue("ale_settings")

        import python_solvers_wrapper_mesh_motion
        self.ale_solver = python_solvers_wrapper_mesh_motion.CreateSolver(self.main_model_part, custom_settings)
        print("::[ALENavierStokesSolverFractionalStep]:: Construction finished")

    def AddVariables(self):
        # add base class variables
        super(ALENavierStokesSolverFractionalStep, self).AddVariables()
        # add mesh motion variables
        self.ale_solver.AddVariables()
        print("::[ALENavierStokesSolverFractionalStep]:: Variables ADDED")

    def AddDofs(self):
        # add base class dofs
        super(ALENavierStokesSolverFractionalStep, self).AddDofs()
        # add mesh motion dofs
        self.ale_solver.AddDofs()
        print("::[ALENavierStokesSolverFractionalStep]:: DOFs ADDED")

    def Initialize(self):
        super(ALENavierStokesSolverFractionalStep, self).Initialize()
        self.ale_solver.Initialize()
        print("::[ALENavierStokesSolverFractionalStep]:: Finished initialization.")

    def GetFluidSolver(self):
        return super(ALENavierStokesSolverFractionalStep, self)

    def GetMeshMotionSolver(self):
        return self.ale_solver

    def Solve(self):
        self.GetMeshMotionSolver().Solve()
        self.GetFluidSolver().Solve()

    def MoveMesh(self):
        self.GetMeshMotionSolver().MoveMesh()
