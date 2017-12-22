from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
CheckForPreviousImport()


def AddVariables(model_part):

    model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(MESH_REACTION)
    model_part.AddNodalSolutionStepVariable(MESH_RHS)

#     print("Mesh solver variables added correctly.")
#
#
def AddDofs(model_part):

    for node in model_part.Nodes:
        node.AddDof(MESH_DISPLACEMENT_X, MESH_REACTION_X)
        node.AddDof(MESH_DISPLACEMENT_Y, MESH_REACTION_Y)
        node.AddDof(MESH_DISPLACEMENT_Z, MESH_REACTION_Z)

import SurfaceTension_monolithic_solver

class MeshSolver:

    def __init__(self, model_part, domain_size, reform_dofs_each_step):

        # Assign parameters
        self.time_order = 2
        self.model_part = model_part
        self.domain_size = domain_size
        self.reform_dofs_each_step = reform_dofs_each_step

        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part, number_of_avg_elems, number_of_avg_nodes)

        # definition of the solvers
        pILUPrecond = ILU0Preconditioner()
        self.linear_solver = BICGSTABSolver(1e-9, 300)
        # self.linear_solver =  DeflatedCGSolver(1e-6, 3000, True,1000)
        # self.linear_solver = ScalingSolver( DeflatedCGSolver(1e-6, 3000, 1000) , True)
        # self.linear_solver = ScalingSolver( DeflatedCGSolver(1e-6, 3000, True,1000) , True)

    def Initialize(self):
        (self.neighbour_search).Execute()
        
                # Initialize mesh model part
        #for node in self.model_part.Nodes:
            #zero = Vector(3)
            #zero[0] = 0.0
            #zero[1] = 0.0
            #zero[2] = 0.0
            #node.SetSolutionStepValue(MESH_REACTION,0,zero)
            #node.SetSolutionStepValue(MESH_DISPLACEMENT,0,zero)
            #node.SetSolutionStepValue(MESH_RHS,0,zero)

        self.solver = ALEApplication.StructuralMeshMovingStrategy(self.model_part, self.linear_solver, self.domain_size, self.time_order, self.reform_dofs_each_step)
        (self.solver).SetEchoLevel(0)
        print("finished moving the mesh")

    def Solve(self):
        if(self.reform_dofs_each_step):
            (self.neighbour_search).Execute()

        #(self.solver).Solve()
        (self.solver).MoveMesh()

    def MoveMesh(self):
        (self.solver).MoveMesh()
        
    #def UpdateReferenceMesh(self):
        #(self.solver).UpdateReferenceMesh()