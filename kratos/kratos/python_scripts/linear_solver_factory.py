from KratosMultiphysics import *

#######################################################################################
#######################################################################################
#######################################################################################
def ConstructPreconditioner( configuration ):
    if hasattr(configuration, 'preconditioner_type'):
      preconditioner_type = configuration.preconditioner_type
      if(preconditioner_type == "None"):
	  return None
      else:
	  if(preconditioner_type == "DiagonalPreconditioner"):
		return DiagonalPreconditioner()
	  elif(preconditioner_type == "ILU0Preconditioner"):
		return ILU0Preconditioner()
	  else:
		print "Preconditioner type not found. Returning None"
		return None
    else:
      return None
      
      


#######################################################################################
#######################################################################################
#######################################################################################
def ConstructSolver( configuration ):
    solver_type = configuration.solver_type
    
    scaling = False
    if hasattr(configuration, 'scaling'):
	scaling = configuration.scaling
    
    linear_solver = None
        
    #######################################################################################
    if(solver_type == "Conjugate gradient"): 
        precond = ConstructPreconditioner( configuration )
        max_it = configuration.max_iteration
	tol    = configuration.tolerance
        if(precond == None):
            linear_solver = CGSolver(tol, max_it)
        else:
            linear_solver = CGSolver(tol, max_it, precond)
    #######################################################################################
    elif(solver_type == "BiConjugate gradient stabilized"):  
        precond = ConstructPreconditioner( configuration )
        max_it = configuration.max_iteration
	tol    = configuration.tolerance
        if(precond == None):
          linear_solver = BICGSTABSolver(tol, max_it)
        else:
          linear_solver = BICGSTABSolver(tol, max_it, precond)
    #######################################################################################
    elif(solver_type == "GMRES"):     
        import KratosMultiphysics.ExternalSolversApplication
        precond = ConstructPreconditioner( configuration )
        max_it = configuration.max_iteration
	tol    = configuration.tolerance
        if(precond == None):
	    linear_solver = KratosMultiphysics.ExternalSolversApplication.GMRESSolver(tol, max_it)
        else:
	    linear_solver = KratosMultiphysics.ExternalSolversApplication.GMRESSolver(tol, max_it, precond)
    #######################################################################################
    elif(solver_type == "Deflated Conjugate gradient"):     
        max_it = configuration.max_iteration
	tol    = configuration.tolerance
	
        assume_constant_structure = False
        if hasattr(configuration, 'assume_constant_structure'):
	    assume_constant_structure = configuration.assume_constant_structure
	  
        max_reduced_size = 1000
        if hasattr(configuration, 'assume_constant_structure'):
	    max_reduced_size = configuration.max_reduced_size
	    
        linear_solver =  DeflatedCGSolver(tol,max_it,assume_constant_structure,max_reduced_size) 
    #######################################################################################
    elif(solver_type == "GMRES-UP Block"):  
	velocity_linear_solver = ConstructSolver( configuration.velocity_block_configuration )
	pressure_linear_solver = ConstructSolver( configuration.pressure_block_configuration )
	m = configuration.gmres_krylov_space_dimension
	max_it = configuration.max_iteration
	tol    = configuration.tolerance
	linear_solver = MixedUPLinearSolver(velocity_linear_solver,pressure_linear_solver,tol,max_it,m)
    #######################################################################################
    elif(solver_type == "Skyline LU factorization"): 
        linear_solver = SkylineLUFactorizationSolver()
    #######################################################################################
    elif(solver_type == "Super LU"):  
	import KratosMultiphysics.ExternalSolversApplication
        linear_solver = KratosMultiphysics.ExternalSolversApplication.SuperLUSolver()
    #######################################################################################
    elif(solver_type == "SuperLUIterativeSolver"):  
	import KratosMultiphysics.ExternalSolversApplication
	tol = configuration.tolerance
	max_it = configuration.max_iteration
	restart = configuration.gmres_krylov_space_dimension
	DropTol = configuration.DropTol
	FillTol = configuration.FillTol
	ilu_level_of_fill = configuration.ilu_level_of_fill
	linear_solver = KratosMultiphysics.ExternalSolversApplication.SuperLUIterativeSolver(tol,max_it,restart,DropTol,FillTol,ilu_level_of_fill)
    #######################################################################################
    elif(solver_type == "PastixDirect"):  
	import KratosMultiphysics.ExternalSolversApplication
	is_symmetric = False
	if hasattr(configuration, 'is_symmetric'):
	  configuration.is_symmetric
	verbosity = 0
	if hasattr(configuration, 'verbosity'):
	    verbosity = configuration.verbosity
        linear_solver = KratosMultiphysics.ExternalSolversApplication.PastixSolver(verbosity,is_symmetric)
    #######################################################################################
    elif(solver_type == "PastixIterative"):  
	import KratosMultiphysics.ExternalSolversApplication
	tol = configuration.tolerance
	max_it = configuration.max_iteration
	restart = configuration.gmres_krylov_space_dimension
	ilu_level_of_fill = configuration.ilu_level_of_fill
	is_symmetric = configuration.is_symmetric
	verbosity = 0
	if hasattr(configuration, 'verbosity'):
	    verbosity = configuration.verbosity
        linear_solver = KratosMultiphysics.ExternalSolversApplication.PastixSolver(tol,restart,ilu_level_of_fill,verbosity,is_symmetric)
    #######################################################################################
    elif(solver_type == "AMGCL"):  
	import KratosMultiphysics.ExternalSolversApplication
	if hasattr(configuration, 'preconditioner_type'):
	  if(configuration.preconditioner_type != "None"):
	    print "WARNING: preconditioner specified in preconditioner_type will not be used as it is not compatible with the AMGCL solver"
	    
        max_it = configuration.max_iteration
	tol    = configuration.tolerance
	
	verbosity = 0
	if hasattr(configuration, 'verbosity'):
	    verbosity = configuration.verbosity
	    
	smoother_type  	= configuration.smoother_type #options are DAMPED_JACOBI, ILU0, SPAI
		
	if(smoother_type == "ILU0"): 
	    amgcl_smoother = KratosMultiphysics.ExternalSolversApplication.AMGCLSmoother.ILU0
	elif(smoother_type == "DAMPED_JACOBI"): 
	    amgcl_smoother = KratosMultiphysics.ExternalSolversApplication.AMGCLSmoother.DAMPED_JACOBI
	elif(smoother_type == "SPAI0"): 
	    amgcl_smoother = KratosMultiphysics.ExternalSolversApplication.AMGCLSmoother.SPAI0
	else:
	    print "ERROR: smoother_type shall be one of ILU0, DAMPED_JACOBI, SPAI0"
	    return None

	krylov_type    	= configuration.krylov_type #options are GMRES, BICGSTAB, CG
	if(krylov_type == "GMRES"): 
	    amgcl_krylov_type = KratosMultiphysics.ExternalSolversApplication.AMGCLIterativeSolverType.GMRES
	elif(krylov_type == "BICGSTAB"): 
	    amgcl_krylov_type = KratosMultiphysics.ExternalSolversApplication.AMGCLIterativeSolverType.BICGSTAB
	elif(krylov_type == "CG"): 
	    amgcl_krylov_type = KratosMultiphysics.ExternalSolversApplication.AMGCLIterativeSolverType.CG
	else:
	    print "ERROR: krylov_type shall be one of GMRES, BICGSTAB, CG"
	    return None

	
	if hasattr(configuration, 'gmres_krylov_space_dimension'):
	    m = configuration.gmres_krylov_space_dimension
	else:
	    m = max_it
        linear_solver = KratosMultiphysics.ExternalSolversApplication.AMGCLSolver(amgcl_smoother,amgcl_krylov_type,tol,max_it,verbosity,m) 
    #######################################################################################
    elif (solver_type == "Parallel MKL Pardiso"):  
	import KratosMultiphysics.MKLSolversApplication
        linear_solver = KratosMultiphysics.MKLSolversApplication.ParallelMKLPardisoSolver()
    else:
	print "*****************************************************************"
	print "Inexisting solver type. Possibilities are:"
	print "Conjugate gradient"
	print "BiConjugate gradient stabilized"
	print "GMRES"
	print "Deflated Conjugate gradient"
	print "AMGCL"
	print "GMRES-UP Block"
	print "Skyline LU factorization"
	print "Super LU (requires ExternalSolversApplication)"
	print "SuperLUIterativeSolver (requires ExternalSolversApplication)"
	print "PastixDirect (requires ExternalSolversApplication + shall be habilitated at compilation time)"
	print "PastixIterative (requires ExternalSolversApplication + shall be habilitated at compilation time)"
	print "Parallel MKL Pardiso (requires MKLSolversApplication)"
	print "*****************************************************************"
	err
    #else:
	#except LogicError:
    
    if(scaling == False):
        return linear_solver
    else:
        return ScalingSolver(linear_solver,True)

