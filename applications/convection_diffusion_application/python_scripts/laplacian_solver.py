from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

def AddVariables(model_part, config=None):
    model_part.AddNodalSolutionStepVariable(TEMPERATURE)



def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(TEMPERATURE)

class LaplacianSolver:
    def __init__(self, model_part, domain_size):
        self.echo_level = 2
        self.model_part = model_part
        self.domain_size = domain_size
        self.linear_solver = None
        self.compute_reactions = False
        self.ReformDofSetAtEachStep = False
        self.MoveMeshFlag = False
        self.calculate_solution_norm = False

    def Initialize(self):
        time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        self.solver = ResidualBasedLinearStrategy(
            self.model_part, time_scheme, self.linear_solver,
            self.compute_reactions, self.ReformDofSetAtEachStep, self.calculate_solution_norm, self.MoveMeshFlag)

        (self.solver).SetEchoLevel(self.echo_level)
        self.solver.Check()



    def Solve(self):
        (self.solver).Solve()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
    def Clear(self):
        (self.solver).Clear()

def CreateSolver(model_part, config):
    solver = LaplacianSolver(model_part, config.domain_size)


    import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
    if(hasattr(config, "linear_solver_config")):
        solver.linear_solver = linear_solver_factory.ConstructSolver(
            config.linear_solver_config)

    return solver
