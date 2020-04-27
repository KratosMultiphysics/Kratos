from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.EigenSolversApplication as KratosEigen
from KratosMultiphysics import python_linear_solver_factory

def ConstructSolver(configuration):
    if not isinstance(configuration, KM.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_type = configuration["solver_type"].GetString()
    if "Application." in solver_type: # the module in which the solver is implemented was specified
        split_name = solver_type.split(".")
        solver_type = split_name[1]

    if KratosEigen.DenseLinearSolverFactory().Has(solver_type):
        KM.Logger.PrintInfo("Eigen-Linear-Solver-Factory", "Constructing a dense linear-solver")
        return KratosEigen.DenseLinearSolverFactory().Create(configuration)
    elif KratosEigen.ComplexDenseLinearSolverFactory().Has(solver_type):
        KM.Logger.PrintInfo("Eigen-Linear-Solver-Factory", "Constructing a complex dense linear-solver")
        return KratosEigen.ComplexDenseLinearSolverFactory().Create(configuration)
    else:
        # fall back to the standard linear solver factory
        return python_linear_solver_factory.ConstructSolver(configuration)
