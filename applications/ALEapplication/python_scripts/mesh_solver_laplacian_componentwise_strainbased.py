from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)


def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(DISPLACEMENT_X)
        node.AddDof(DISPLACEMENT_Y)
        node.AddDof(DISPLACEMENT_Z)

    print("variables for the mesh solver added correctly")


class MeshSolver:

    def __init__(self, model_part, domain_size, reform_dof_at_every_step):

        # Assign parameters
        self.model_part = model_part
        self.domain_size = domain_size
        self.reform_dof_at_every_step = reform_dof_at_every_step
        self.time_order = 2

        # Neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part, number_of_avg_elems, number_of_avg_nodes)

        # definition of the solvers
        pILUPrecond = ILU0Preconditioner()
        self.linear_solver = DeflatedCGSolver(1e-6, 3000, True, 1000)
        # self.linear_solver =  BICGSTABSolver(1e-3, 300,pILUPrecond)
        self.linear_solver = ScalingSolver(DeflatedCGSolver(1e-6, 3000, True, 1000), True)

    def Initialize(self):
        (self.neighbour_search).Execute()

        self.solver = LaplacianComponentwiseMeshMovingStrategyStrainbased(self.model_part, self.linear_solver, self.domain_size, self.time_order, self.reform_dof_at_every_step)
        (self.solver).SetEchoLevel(0)
        print("finished imoving the mesh")

    def Solve(self):
        if(self.reform_dof_at_every_step):
            (self.neighbour_search).Execute()

        (self.solver).Solve()

    def MoveNodes(self):
        (self.solver).MoveNodes()
