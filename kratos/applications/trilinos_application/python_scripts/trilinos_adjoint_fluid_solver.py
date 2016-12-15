from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.mpi import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(ADJOINT_VELOCITY)
    model_part.AddNodalSolutionStepVariable(ADJOINT_PRESSURE)
    model_part.AddNodalSolutionStepVariable(PRIMAL_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRIMAL_PRESSURE)
    model_part.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)
    model_part.AddNodalSolutionStepVariable(NORMAL_SENSITIVITY)
    mpi.world.barrier()
    if mpi.rank == 0:
        print("Variables for the trilinos adjoint fluid solver added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(ADJOINT_VELOCITY_X)
        node.AddDof(ADJOINT_VELOCITY_Y)
        node.AddDof(ADJOINT_VELOCITY_Z)
        node.AddDof(ADJOINT_PRESSURE)
    mpi.world.barrier()
    if mpi.rank == 0:
        print("Dofs for the trilinos adjoint fluid solver added correctly")


class TrilinosAdjointFluidSolver:
    
    def __init__(self, model_part, dimension):

        self.model_part = model_part
        self.dimension = dimension
        self.linear_solver_tol = 1e-6
        self.linear_solver_max_it = 1000

        aztec_parameters = ParameterList()
        aztec_parameters.set("AZ_solver", "AZ_gmres")
        aztec_parameters.set("AZ_output", "AZ_none")
        MLList = ParameterList()
        default_settings = EpetraDefaultSetter()
        default_settings.SetDefaults(MLList, "NSSA")
        MLList.set("ML output", 0)
        MLList.set("max levels", 3)
        MLList.set("aggregation: type", "Uncoupled")
        self.linear_solver = MultiLevelSolver(aztec_parameters, MLList, self.linear_solver_tol, self.linear_solver_max_it)

        self.comm = CreateCommunicator()
        
        mpi.world.barrier()
        if mpi.rank == 0:
            print("Construction adjoint fluid solver finished")

    def Initialize(self):
        self.solver = TrilinosAdjointFluidStrategy(self.comm,
                                                   self.model_part,
                                                   self.linear_solver,
                                                   self.dimension)
        mpi.world.barrier()
        if mpi.rank == 0:
            print ("Initialization adjoint fluid solver finished")
    
    def Solve(self):
        (self.solver).Solve()

    def SetDragForceDirection(self, direction):
        (self.solver).SetDragForceDirection(direction)

    def ComputeSensitivity(self):
        (self.solver).ComputeSensitivity()
