from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7



def ConstructPreconditioner(configuration):
    import KratosMultiphysics
    if(configuration["preconditioner_type"].GetString() == "None"):
        return KratosMultiphysics.Preconditioner()
    elif(configuration["preconditioner_type"].GetString() == "DiagonalPreconditioner"):
        return KratosMultiphysics.DiagonalPreconditioner()
    elif(configuration["preconditioner_type"].GetString() == "ILU0Preconditioner"):
        return KratosMultiphysics.ILU0Preconditioner()
    elif(configuration["preconditioner_type"].GetString() == "ILUPreconditioner"):
        return KratosMultiphysics.ILUPreconditioner() 
    else:
        raise Exception("preconditioner_type specified is not correct. Current choice is :" + configuration["preconditioner_type"].GetString())

def ConstructSolver(configuration):
    import KratosMultiphysics
    
    if(type(configuration) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")
    
    solver_type = configuration["solver_type"].GetString()
    
    if configuration.Has("scaling"):
        scaling = configuration["scaling"].GetBool()
    else:
        scaling = False
        
    if(solver_type == "CGSolver"):
        linear_solver = KratosMultiphysics.CGSolver( configuration, ConstructPreconditioner(configuration) )
    elif(solver_type == "BICGSTABSolver"):
        linear_solver = KratosMultiphysics.BICGSTABSolver( configuration, ConstructPreconditioner(configuration) )
    elif(solver_type == "TFQMRSolver"):
        linear_solver = KratosMultiphysics.TFQMRSolver( configuration, ConstructPreconditioner(configuration) )
    elif(solver_type == "DeflatedCGSolver"):
        linear_solver = KratosMultiphysics.DeflatedCGSolver( configuration )
    elif(solver_type == "MixedUPLinearSolver"):
        linear_solver = KratosMultiphysics.MixedUPLinearSolver(configuration)
    elif(solver_type == "SkylineLUFactorizationSolver"):
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver(configuration)
    elif(solver_type == "complex_skyline_lu_solver"):
        linear_solver = KratosMultiphysics.ComplexSkylineLUSolver(configuration)

    ################################## following solvers need importing the ExternalSolversApplication
    elif(solver_type == "GMRESSolver"):
        import KratosMultiphysics.ExternalSolversApplication
        linear_solver = KratosMultiphysics.ExternalSolversApplication.GMRESSolver( configuration, ConstructPreconditioner(configuration) )
    elif(solver_type == "SuperLUSolver" or solver_type=="Super LU" or solver_type=="Super_LU"):
        import KratosMultiphysics.ExternalSolversApplication
        linear_solver = KratosMultiphysics.ExternalSolversApplication.SuperLUSolver(configuration)
    elif(solver_type == "SuperLUIterativeSolver"):
        import KratosMultiphysics.ExternalSolversApplication
        linear_solver = KratosMultiphysics.ExternalSolversApplication.SuperLUIterativeSolver(configuration)
    elif(solver_type == "PastixSolver"):
        import KratosMultiphysics.ExternalSolversApplication
        linear_solver = KratosMultiphysics.ExternalSolversApplication.PastixSolver(configuration)
    elif(solver_type == "AMGCL"):
        linear_solver = KratosMultiphysics.AMGCLSolver(configuration)
    elif(solver_type == "AMGCL_NS_Solver"):
        linear_solver = KratosMultiphysics.AMGCL_NS_Solver(configuration)
    elif(solver_type == "complex_pastix_solver"):
        import KratosMultiphysics.ExternalSolversApplication
        linear_solver = KratosMultiphysics.ExternalSolversApplication.PastixComplexSolver(configuration)
 
    ################################## following solvers need importing the MKLSolversApplication
    elif (solver_type == "ParallelMKLPardisoSolver"):
        import KratosMultiphysics.MKLSolversApplication
        linear_solver = KratosMultiphysics.MKLSolversApplication.ParallelMKLPardisoSolver(configuration)

    ###################################### FAILED TO FIND solver_type
    else:
        raise Exception("solver type not found. Asking for :" + solver_type)


    ###### here decide if a prescaling is to be applied 
    if(scaling == False):
        return linear_solver
    else:
        return KratosMultiphysics.ScalingSolver(linear_solver, True)
