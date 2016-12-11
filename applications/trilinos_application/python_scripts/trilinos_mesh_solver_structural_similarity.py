from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.TrilinosApplication import *
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
        node.AddDof(MESH_DISPLACEMENT_X, MESH_REACTION_X)
        node.AddDof(MESH_DISPLACEMENT_Y, MESH_REACTION_Y)
        node.AddDof(MESH_DISPLACEMENT_Z, MESH_REACTION_Z)
    mpi.world.barrier()
    if mpi.rank == 0:
        print("dofs for the mesh solver added correctly")


class TrilinosMeshSolverStructuralSimilarity:

    def __init__(self, model_part, domain_size, reform_dof_at_every_step=False):

        # Assigning parameters to be used
        self.time_order = 2
        self.model_part = model_part
        self.domain_size = domain_size
        self.reform_dof_at_every_step = reform_dof_at_every_step

        # Create communicator
        self.Comm = CreateCommunicator()

        # Define solver
        import trilinos_linear_elastic_ml_solver
        nit_max = 10000
        linear_tol = 1e-5
        self.linear_solver = trilinos_linear_elastic_ml_solver.MultilevelLinearSolver(linear_tol, nit_max)

    def Initialize(self):
        # Initialize mesh model part
        for node in self.model_part.Nodes:
            zero = Vector(3)
            zero[0] = 0.0
            zero[1] = 0.0
            zero[2] = 0.0
            node.SetSolutionStepValue(MESH_REACTION,0,zero)
            node.SetSolutionStepValue(MESH_DISPLACEMENT,0,zero)
            node.SetSolutionStepValue(MESH_RHS,0,zero)

        self.solver = TrilinosStructuralMeshMovingStrategy(self.Comm, self.model_part, self.linear_solver, self.time_order, self.reform_dof_at_every_step)
        (self.solver).SetEchoLevel(0)
        print("finished moving the mesh")

    def Solve(self):
        (self.solver).Solve()

    def MoveNodes(self):
        (self.solver).MoveNodes()
