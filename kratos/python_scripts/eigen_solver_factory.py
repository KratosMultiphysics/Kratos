from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

def ConstructSolver(settings):
    import KratosMultiphysics
    
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Input is expected to be provided as a Kratos Parameters object")
    
    solver_type = settings["solver_type"].GetString()
    
    if solver_type == "eigen_sparse_eigensystem":
        import KratosMultiphysics.EigenSolversApplication
        eigen_solver = KratosMultiphysics.EigenSolversApplication.SparseEigensystemSolver(settings)
        return eigen_solver
    
    import new_linear_solver_factory
    linear_solver = new_linear_solver_factory.ConstructSolver(settings["linear_solver_settings"])
        
    if(solver_type == "power_iteration_eigenvalue_solver"):
        eigen_solver = KratosMultiphysics.PowerIterationEigenvalueSolver( settings, linear_solver)
    elif(solver_type == "power_iteration_highest_eigenvalue_solver"):
        eigen_solver = KratosMultiphysics.PowerIterationHighestEigenvalueSolver( settings, linear_solver)
    elif(solver_type == "rayleigh_quotient_iteration_eigenvalue_solver"):
        eigen_solver = KratosMultiphysics.RayleighQuotientIterationEigenvalueSolver( settings, linear_solver)
    elif(solver_type == "FEAST" or solver_type == "feast"):
        import KratosMultiphysics.ExternalSolversApplication
        eigen_solver = KratosMultiphysics.ExternalSolversApplication.FEASTSolver(settings, linear_solver)
    else:
        raise Exception("Solver type not found. Asking for :" + solver_type)
    
    return eigen_solver
