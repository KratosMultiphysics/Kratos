from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

def ConstructSolver(configuration):
    import KratosMultiphysics

    depr_msg  = '"kratos/python_scripts/new_linear_solver_factory.py" is deprecated and will be removed\n'
    depr_msg += 'Please use "kratos/python_scripts/python_linear_solver_factory.py" instead!'
    KratosMultiphysics.Logger.PrintWarning('DEPRECATION-WARNING', depr_msg)

    if(type(configuration) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    if KratosMultiphysics.ComplexLinearSolverFactory().Has(configuration["solver_type"].GetString()):
        return KratosMultiphysics.ComplexLinearSolverFactory().Create(configuration)
    else:
        return KratosMultiphysics.LinearSolverFactory().Create(configuration)
