from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MeshMovingApplication", "FluidDynamicsApplication","TrilinosApplication")

# Import applications
import KratosMultiphysics.MeshMovingApplication as MeshMovingApplication
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

# Other imports
import trilinos_navier_stokes_solver_vmsmonolithic
#import trilinos_mesh_solver_base


def CreateSolver(model_part, custom_settings):
    return ALETrilinosNavierStokesSolverVMSMonolithic(model_part, custom_settings)


class ALETrilinosNavierStokesSolverVMSMonolithic(trilinos_navier_stokes_solver_vmsmonolithic.TrilinosNavierStokesSolverMonolithic):
    def __init__(self, model_part, custom_settings):
        # remove the ale_settings so we can use the navier stokes constructor
        navier_stokes_settings = custom_settings.Clone()
        navier_stokes_settings.RemoveValue("ale_settings")
        navier_stokes_settings["solver_type"].SetString("Monolithic")
        super(ALETrilinosNavierStokesSolverVMSMonolithic, self).__init__(model_part, navier_stokes_settings)
        # create mesh motion solver
        custom_settings.AddEmptyValue("problem_data")
        custom_settings["problem_data"].AddEmptyValue("parallel_type")
        custom_settings["problem_data"]["parallel_type"].SetString("MPI")
        custom_settings.AddValue("solver_settings", custom_settings["ale_settings"])
        custom_settings.RemoveValue("ale_settings")

        import python_solvers_wrapper_mesh_motion
        self.ale_solver = python_solvers_wrapper_mesh_motion.CreateSolver(self.main_model_part, custom_settings)
        if KratosMPI.mpi.rank == 0:
            print("::[ALETrilinosNavierStokesSolverVMSMonolithic]:: Construction finished")

    def AddVariables(self):
        # add base class variables
        super(ALETrilinosNavierStokesSolverVMSMonolithic, self).AddVariables()
        self.ale_solver.AddVariables()
        if KratosMPI.mpi.rank == 0:
            print("::[ALETrilinosNavierStokesSolverVMSMonolithic]:: Variables ADDED.")

    def AddDofs(self):
        super(ALETrilinosNavierStokesSolverVMSMonolithic, self).AddDofs()
        self.ale_solver.AddDofs()
        if KratosMPI.mpi.rank == 0:
            print("::[ALETrilinosNavierStokesSolverVMSMonolithic]:: DOFs ADDED.")

    def Initialize(self):
        super(ALETrilinosNavierStokesSolverVMSMonolithic, self).Initialize()
        self.ale_solver.Initialize()
        if KratosMPI.mpi.rank == 0:
            print("::[ALETrilinosNavierStokesSolverVMSMonolithic]:: Finished initialization.")

    def GetFluidSolver(self):
        return super(ALETrilinosNavierStokesSolverVMSMonolithic, self)

    def GetMeshMotionSolver(self):
        return self.ale_solver

    def Solve(self):
        self.GetMeshMotionSolver().Solve()
        self.GetFluidSolver.Solve()

    def MoveMesh(self):
        self.GetMeshMotionSolver().MoveMesh()
