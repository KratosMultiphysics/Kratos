from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MeshMovingApplication")

# Import applications
#import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving

# Other imports
#import os
import mesh_solver_base



def CreateSolver(model_part, custom_settings):
    return MeshSolverBallvertex(model_part,custom_settings)


class MeshSolverBallvertex(mesh_solver_base.MeshSolverBase):
    def __init__(self, model_part, custom_settings):
        super(MeshSolverBallvertex, self).__init__(model_part, custom_settings)
        print("::[MeshSolverBallvertex]:: Construction finished")

    #### Public user interface functions ####

    def AddVariables(self):
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        print("::[MeshSolverBallvertex]:: Variables ADDED.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z)
        print("::[MeshSolverBallvertex]:: DOFs ADDED.")

    def Initialize(self):
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        neighbour_search = FindNodalNeighboursProcess(self.model_part,
                                                      number_of_avg_elems,
                                                      number_of_avg_nodes)
        neighbour_search.Execute()
        solver = self.get_mesh_motion_solving_strategy()
        if self.settings["reform_dofs_each_step"].GetBool() == False:
            solver.ConstructSystem(self.model_part)
        print("::[MeshSolverBallvertex]:: Finished initialization.")

    def Solve(self):
        solver = self.get_mesh_motion_solving_strategy()
        linear_solver = self.get_linear_solver()
        if self.settings["reform_dofs_each_step"].GetBool() == True:
            (self.neighbour_search).Execute()
            solver.ConstructSystem(self.model_part)
            solver.BuildAndSolveSystem(self.model_part, linear_solver)
            solver.ClearSystem()
        else:
            solver.BuildAndSolveSystem(model_part, linear_solver)
        # Move mesh and calculate mesh velocity.
        time_order = self.settings("time_order").GetInt()
        MeshMovingApplication.MoveMeshUtilities().BDF_MoveMesh(time_order, self.model_part)

    def MoveMesh(self):
        time_order = self.settings("time_order").GetInt()
        MeshMovingApplication.MoveMeshUtilities().BDF_MoveMesh(time_order, self.model_part)

    #### Private functions ####

    def _create_linear_solver(self):
        linear_solver = ScalingSolver(DeflatedCGSolver(1e-6, 3000, True, 1000), True)
        return linear_solver

    def _create_mesh_motion_solving_strategy(self):
        domain_size = self.model_part.ProcessInfo[DOMAIN_SIZE]
        if(domain_size == 2):
            solver = BallVertexMeshMoving2D()
        else:
            solver = BallVertexMeshMoving3D()
        return solver
