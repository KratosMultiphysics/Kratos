from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

def ConstructSolver(configuration):
    import KratosMultiphysics

    if(type(configuration) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    if KratosMultiphysics.ComplexLinearSolverFactory().Has(configuration["solver_type"].GetString()):
        return KratosMultiphysics.ComplexLinearSolverFactory().CreateSolver(configuration)
    else:
        return KratosMultiphysics.LinearSolverFactory().CreateSolver(configuration)
