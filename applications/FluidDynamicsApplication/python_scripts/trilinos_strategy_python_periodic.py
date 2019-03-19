from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *

# import base class
from trilinos_strategy_python import SolvingStrategyPython


class SolvingStrategyPeriodic(SolvingStrategyPython):

    def __init__(
        self, domain_size, model_part, time_scheme, linear_solver, convergence_criteria,
            CalculateReactionsFlag, ReformDofSetAtEachStep, MoveMeshFlag, Comm, guess_row_size, periodic_var):
        # add my own builder and solver
        self.builder_and_solver = TrilinosBlockBuilderAndSolverPeriodic(
            Comm, guess_row_size, linear_solver, periodic_var)
        # initialize skiping builder and solver
        SolvingStrategyPython.__init__(
            self,
            "NONE",
            model_part,
            time_scheme,
            linear_solver,
            convergence_criteria,
            CalculateReactionsFlag,
            ReformDofSetAtEachStep,
            MoveMeshFlag,
            Comm,
            guess_row_size)
