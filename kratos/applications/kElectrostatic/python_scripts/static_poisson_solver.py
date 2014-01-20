from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from Kratos import *
from KratosR1ElectrostaticApplication import *


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(ELECTRICAL_PERMITTIVITY)
    model_part.AddNodalSolutionStepVariable(ELECTRIC_FIELD)
    model_part.AddNodalSolutionStepVariable(ELECTRIC_DISPLACEMENT_FIELD)
    model_part.AddNodalSolutionStepVariable(ELECTROSTATIC_POTENTIAL)
    model_part.AddNodalSolutionStepVariable(ELECTROSTATIC_POINT_CHARGE)
    model_part.AddNodalSolutionStepVariable(INFINIT_COEFFICIENT)


def AddDofs(model_part):
    for node in model_part.Nodes:

        # adding dofs
        node.AddDof(ELECTROSTATIC_POTENTIAL)

    print("variables for the Poisson solver added correctly")


class StaticPoissonSolver:
    #

    def __init__(self, model_part, domain_size):

        self.model_part = model_part
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        # definition of the solvers
        # self.poisson_linear_solver =  SkylineLUFactorizationSolver()

        pDiagPrecond = DiagonalPreconditioner()
        self.poisson_linear_solver = BICGSTABSolver(1e-6, 5000, pDiagPrecond)

        # definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(1e-6, 1e-9)

    #
    def Initialize(self):
        # creating the solution strategy
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = False
        MoveMeshFlag = False
        import strategy_python
        self.solver = strategy_python.SolvingStrategyPython(self.model_part, self.time_scheme, self.poisson_linear_solver, self.conv_criteria, CalculateReactionFlag, ReformDofSetAtEachStep, MoveMeshFlag)

    #
    def Solve(self):
        (self.solver).Solve()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)
