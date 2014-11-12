from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.TrilinosApplication import *
#from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.mpi import *
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
    mpi.world.barrier()
    if mpi.rank == 0:
        print("variables for the mesh solver added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
        node.AddDof(DISPLACEMENT_Z);
    print("variables for the mesh solver added correctly")

    mpi.world.barrier()
    if mpi.rank == 0:
        print("dofs for the mesh solver added correctly")


class TrilinosMeshSolverStructuralSimilarityNonlin:

    def __init__(self, model_part, domain_size, reform_dof_at_every_step, NonLinearTol, MaxIter):

        # Assignation of parameters to be used
        self.time_order = 2
        self.model_part = model_part
        self.domain_size = domain_size
        self.reform_dof_at_every_step = reform_dof_at_every_step

        # for non-linear
        self.NonLinearTol = NonLinearTol
        self.MaxIter = MaxIter

                # Add dofs
                # AddDofs(model_part)

        # Create communicator
        self.Comm = CreateCommunicator()

        # Define solver
        import MonolithicMultiLevelSolver
        pressure_nit_max = 1000
        pressure_linear_tol = 1e-6
        self.linear_solver = MonolithicMultiLevelSolver.LinearSolver(pressure_linear_tol, pressure_nit_max)

    def Initialize(self):
        #(self.neighbour_search).Execute()

        self.solver = TrilinosStructuralMeshMovingStrategyNonlin(
            self.Comm,
            self.model_part,
            self.linear_solver,
            self.domain_size,
            self.time_order,
            self.reform_dof_at_every_step,
            self.NonLinearTol,
            self.MaxIter)
        (self.solver).SetEchoLevel(0)
        print("finished moving the mesh")

    def Solve(self):
        if(self.reform_dof_at_every_step):
            (self.neighbour_search).Execute()

        (self.solver).Solve()

    def MoveNodes(self):
        (self.solver).MoveNodes()
