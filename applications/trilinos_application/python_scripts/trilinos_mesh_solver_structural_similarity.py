from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.mpi import *
CheckForPreviousImport()

def CreateSolver(model_part, custom_settings):
    return TrilinosMeshSolverStructuralSimilarity(model_part, custom_settings)

class TrilinosMeshSolverStructuralSimilarity:

    def __init__(self, model_part, custom_settings):

        default_settings = Parameters("""
        {
            "time_order"                 : 2,
            "mesh_reform_dofs_each_step" : false
        }""")

        custom_settings.ValidateAndAssignDefaults(default_settings)

        # set parameters
        self.model_part = model_part
        self.time_order = custom_settings["time_order"].GetInt()
        self.mesh_reform_dofs_each_step = custom_settings["mesh_reform_dofs_each_step"].GetBool()

        # Create communicator
        self.Comm = CreateCommunicator()

        # Define solver
        import trilinos_linear_elastic_ml_solver
        nit_max = 10000
        linear_tol = 1e-5
        self.linear_solver = trilinos_linear_elastic_ml_solver.MultilevelLinearSolver(linear_tol, nit_max)
        if mpi.rank == 0:
            print("Construction of TrilinosMeshSolverStructuralSimilarity finished")

    def AddVariables(self):
        self.model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(MESH_REACTION)
        self.model_part.AddNodalSolutionStepVariable(MESH_RHS)
        self.model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
        mpi.world.barrier()
        if mpi.rank == 0:
            print("variables for the mesh solver added correctly")

    def AddDofs(self):
        for node in self.model_part.Nodes:
            node.AddDof(MESH_DISPLACEMENT_X, MESH_REACTION_X)
            node.AddDof(MESH_DISPLACEMENT_Y, MESH_REACTION_Y)
            node.AddDof(MESH_DISPLACEMENT_Z, MESH_REACTION_Z)
        mpi.world.barrier()
        if mpi.rank == 0:
            print("dofs for the mesh solver added correctly")

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

        self.solver = TrilinosStructuralMeshMovingStrategy(self.Comm, self.model_part, self.linear_solver, self.time_order, self.mesh_reform_dofs_each_step)
        (self.solver).SetEchoLevel(0)

    def Solve(self):
        (self.solver).Solve()

    def MoveNodes(self):
        (self.solver).MoveNodes()
