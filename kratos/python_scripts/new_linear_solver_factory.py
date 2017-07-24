from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7


def ConstructSolver(configuration):
    import KratosMultiphysics
    return KratosMultiphysics.LinearSolverFactoryBase().CreateSolver(configuration)
    
 