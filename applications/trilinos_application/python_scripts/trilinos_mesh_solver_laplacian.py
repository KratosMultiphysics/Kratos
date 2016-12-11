from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.mpi import *
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(MESH_REACTION)
    model_part.AddNodalSolutionStepVariable(MESH_RHS)
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
    mpi.world.barrier()
    if mpi.rank == 0:
        print("variables for the mesh solver added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(MESH_DISPLACEMENT_X, MESH_REACTION_X)
        node.AddDof(MESH_DISPLACEMENT_Y, MESH_REACTION_Y)
        node.AddDof(MESH_DISPLACEMENT_Z, MESH_REACTION_Z)

    mpi.world.barrier()
    if mpi.rank == 0:
        print("dofs for the mesh solver added correctly")


class TrilinosMeshSolverComponentwise:

    def __init__(self, model_part, domain_size, reform_dof_at_every_step):

        self.model_part = model_part
        self.domain_size = domain_size
        self.reform_dof_at_every_step = reform_dof_at_every_step

        #AddDofs(model_part)

        # assignation of parameters to be used
        self.time_order = 1

        # Create communicator
        self.Comm = CreateCommunicator()

        # Define solver
        import PressureMultiLevelSolver
        pressure_nit_max = 1000
        pressure_linear_tol = 1e-6
        self.linear_solver = PressureMultiLevelSolver.MultilevelLinearSolver(pressure_linear_tol, pressure_nit_max)

    def Initialize(self):
        self.solver = TrilinosLaplacianMeshMovingStrategy(self.Comm, self.model_part, self.linear_solver, self.domain_size, self.time_order, self.reform_dof_at_every_step)
        (self.solver).SetEchoLevel(0)
        print("finished moving the mesh")

    def Solve(self):
        (self.solver).Solve()

    def MoveNodes(self):
        (self.solver).MoveNodes()
