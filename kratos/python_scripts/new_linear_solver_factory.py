def ConstructSolver(configuration):
    import KratosMultiphysics

    depr_msg  = '"kratos/python_scripts/new_linear_solver_factory.py" is deprecated and will be removed\n'
    depr_msg += 'Please use "kratos/python_scripts/python_linear_solver_factory.py" instead!'
    KratosMultiphysics.Logger.PrintWarning('DEPRECATION-WARNING', depr_msg)

    from KratosMultiphysics import python_linear_solver_factory
    return python_linear_solver_factory.ConstructSolver(configuration)
