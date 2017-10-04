from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

def ConstructSolver(settings):
    import KratosMultiphysics
    
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Input is expected to be provided as a Kratos Parameters object")
    
    default_settings = KratosMultiphysics.Parameters("""
    {
        "eigen_solver_settings"      : {
            "solver_type"             : "PowerIterationEigenvalueSolver",
            "max_iteration"           : 10000,
            "tolerance"               : 1e-8,
            "required_eigen_number"   : 1,
            "shifting_convergence"    : 0.25,
            "verbosity"               : 0
        },
        "linear_solver_settings"      : {
            "solver_type"             : "SuperLUSolver",
            "max_iteration"           : 500,
            "tolerance"               : 1e-9,
            "scaling"                 : false,
            "verbosity"               : 0
        }
    }
    """
    )
    
    settings.ValidateAndAssignDefaults(default_settings)
    
    import new_linear_solver_factory
    linear_solver = new_linear_solver_factory.ConstructSolver(settings["linear_solver_settings"])
    
    solver_type = settings["eigen_solver_settings"]["solver_type"].GetString()
        
    if(solver_type == "PowerIterationEigenvalueSolver"):
        eigen_solver = KratosMultiphysics.PowerIterationEigenvalueSolver( settings["eigen_solver_settings"], linear_solver)
    elif(solver_type == "PowerIterationHighestEigenvalueSolver"):
        eigen_solver = KratosMultiphysics.PowerIterationHighestEigenvalueSolver( settings["eigen_solver_settings"], linear_solver)
    elif(solver_type == "RayleighQuotientIterationEigenvalueSolver"):
        eigen_solver = KratosMultiphysics.RayleighQuotientIterationEigenvalueSolver( settings["eigen_solver_settings"], linear_solver)

    ################################### Following solvers need importing the ExternalSolversApplication # TODO: Ask Michael Andre to unify interface
    #elif(solver_type == "FEASTSolver"): 
        #import KratosMultiphysics.ExternalSolversApplication
        #eigen_solver = KratosMultiphysics.ExternalSolversApplication.FEASTSolver( settings["eigen_solver_settings"], linear_solver )

    ###################################### FAILED TO FIND solver_type
    else:
        raise Exception("Solver type not found. Asking for :" + solver_type)
    
    return eigen_solver
