from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(ELEMENTSHAPE)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)


def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(DISPLACEMENT_X)
        node.AddDof(DISPLACEMENT_Y)
        node.AddDof(DISPLACEMENT_Z)
    print("variables for the mesh solver added correctly")


class MeshSolverStructuralSimilarityNonlinear:

    def __init__(self, model_part, domain_size, reform_dof_at_every_step, NonLinearTol, MaxIter):

        # set parameters
        self.time_order = 2
        self.model_part = model_part
        self.domain_size = domain_size
        self.reform_dof_at_every_step = reform_dof_at_every_step

        # for non-linear
        self.NonLinearTol = NonLinearTol
        self.MaxIter = MaxIter

        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part, number_of_avg_elems, number_of_avg_nodes)

        # definition of the solvers
        tol = 1e-6
        max_it = 500
        verbosity = 1
        m = 100
        self.linear_solver = AMGCLSolver(AMGCLSmoother.DAMPED_JACOBI, AMGCLIterativeSolverType.BICGSTAB, tol, max_it, verbosity, m)


    def Initialize(self):
        (self.neighbour_search).Execute()

        self.solver = StructuralMeshMovingStrategyNonlinear(self.model_part, self.linear_solver, self.domain_size, self.time_order, self.reform_dof_at_every_step, self.NonLinearTol, self.MaxIter)
        (self.solver).SetEchoLevel(0)
        print("Finished moving the mesh")

    def Solve(self):
        if(self.reform_dof_at_every_step):
            (self.neighbour_search).Execute()

        (self.solver).Solve()

    def MoveNodes(self):
        (self.solver).MoveNodes()
