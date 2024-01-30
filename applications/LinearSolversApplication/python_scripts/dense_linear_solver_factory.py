
import KratosMultiphysics as KM
import KratosMultiphysics.LinearSolversApplication as KratosSolvers
from KratosMultiphysics import python_linear_solver_factory

def ConstructSolver(configuration):
    if not isinstance(configuration, KM.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_type = configuration["solver_type"].GetString()
    if "Application." in solver_type: # the module in which the solver is implemented was specified
        split_name = solver_type.split(".")
        solver_type = split_name[1]

    if KratosSolvers.DenseLinearSolverFactory().Has(solver_type):
        KM.Logger.PrintInfo("Eigen-Linear-Solver-Factory", "Constructing a dense linear-solver")
        return KratosSolvers.DenseLinearSolverFactory().Create(configuration)
    elif KratosSolvers.ComplexDenseLinearSolverFactory().Has(solver_type):
        KM.Logger.PrintInfo("Eigen-Linear-Solver-Factory", "Constructing a complex dense linear-solver")
        return KratosSolvers.ComplexDenseLinearSolverFactory().Create(configuration)
    else:
        # fall back to the standard linear solver factory
        return python_linear_solver_factory.ConstructSolver(configuration)
