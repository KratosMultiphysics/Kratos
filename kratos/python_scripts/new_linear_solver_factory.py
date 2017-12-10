from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7


def ConstructSolver(configuration):
    import KratosMultiphysics
    if KratosMultiphysics.ComplexLinearSolverFactoryBase().Has(configuration["solver_type"].GetString()):
        return KratosMultiphysics.ComplexLinearSolverFactoryBase().CreateSolver(configuration)
    else:
        return KratosMultiphysics.LinearSolverFactoryBase().CreateSolver(configuration)
    
