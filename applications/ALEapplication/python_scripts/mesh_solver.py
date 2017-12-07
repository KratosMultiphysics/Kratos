from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
import KratosMultiphysics.ALEApplication as ALEApplication
CheckForPreviousImport()
import mesh_solver



def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(AUX_MESH_VAR)


def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(AUX_MESH_VAR)
    print("variables for the mesh solver added correctly")


def CreateSolver(model_part, domain_size, reform_dof_at_every_ste):
    return MeshSolver(model_part, domain_size, reform_dof_at_every_ste)


class MeshSolver:

    def __init__(self, model_part, domain_size, reform_dof_at_every_step):

        # Assign parameters
        self.time_order = 2
        self.model_part = model_part
        self.domain_size = 2
        self.reform_dof_at_every_step = reform_dof_at_every_step

        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part, number_of_avg_elems, number_of_avg_nodes)
        
        self.compute_reactions = False


        # definition of the solvers
        pILUPrecond = ILU0Preconditioner()
        self.linear_solver = BICGSTABSolver(1e-9, 300)

        
    def AddVariables(self):
        super(MeshSolver, self).AddVariables()
        self.ale_solver.AddVariables()
        print("::[MeshSolver]:: Variables ADDED.")

    def AddDofs(self):
        super(MeshSolver, self).AddDofs()
        self.ale_solver.AddDofs()
        print("::[MeshSolver]:: DOFs ADDED.")
        

    def Initialize(self):
        (self.neighbour_search).Execute()

        self.solver = ALEApplication.StructuralMeshMovingStrategy(self.model_part, self.linear_solver, self.time_order, self.reform_dof_at_every_step, self.compute_reactions)
        (self.solver).SetEchoLevel(0)
        print("finished moving the mesh")



    ##def GetFluidSolver(self):
        ##return super(MeshSolver, self)

    def GetMeshMotionSolver(self):
        return self.ale_solver

    def Solve(self):
        if(self.reform_dof_at_every_step):
            (self.neighbour_search).Execute()
        self.GetMeshMotionSolver().Solve()

    def MoveMesh(self):
        self.GetMeshMotionSolver().MoveMesh()

