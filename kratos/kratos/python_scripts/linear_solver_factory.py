from KratosMultiphysics import *

#######################################################################################
#######################################################################################
#######################################################################################
def ConstructPreconditioner( configuration ):
    try:
      preconditioner_type = configuration.preconditioner_type
    except LogicError:
      print "missing preconditioner type. Either set the preconditioner type as None or to a given preconditioner"
      
      
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

    
#######################################################################################
#######################################################################################
#######################################################################################
def ConstructSolver( configuration ):
    solver_type = configuration.solver_type
    max_it = configuration.max_iteration
    tol    = configuration.tolerance
    scaling = configuration.scaling
    
    linear_solver = None
        
    #######################################################################################
    if(solver_type == "Conjugate gradient"): 
        precond = ConstructPreconditioner( configuration )
        if(precond == None):
            linear_solver = CGSolver(tol, max_it)
        else:
            linear_solver = CGSolver(tol, max_it, precond)
    #######################################################################################
    elif(solver_type == "BiConjugate gradient stabilized"):  
        precond = ConstructPreconditioner( configuration )
        if(precond == None):
          linear_solver = BICGSTABSolver(tol, max_it)
        else:
          linear_solver = BICGSTABSolver(tol, max_it, precond)
    #######################################################################################
    elif(solver_type == "GMRES"):     
        precond = ConstructPreconditioner( configuration )
        if(precond == None):
	    linear_solver = GMRESSolver(tol, max_it)
        else:
	    linear_solver = GMRESSolver(tol, max_it, precond)
    #######################################################################################
    elif(solver_type == "Mixed UP"):  
	velocity_linear_solver = ConstructSolver( configuration.velocity_linear_solver_configuration )
	pressure_linear_solver = ConstructSolver( configuration.pressure_linear_solver_configuration )
	m = configuration.gmres_krylov_space_dimension
	linear_solver = MixedUPLinearSolver(velocity_linear_solver,pressure_linear_solver,tol,max_it,m)
    #######################################################################################
    elif(solver_type == "Skyline LU factorization"): 
        linear_solver = SkylineLUFactorizationSolver()
    #######################################################################################
    elif(solver_type == "Super LU"):  
	import KratosMultiphysics.ExternalSolversApplication
        linear_solver = KratosMultiphysics.ExternalSolversApplication.SuperLUSolver()
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
	print "Mixed UP"
	print "Skyline LU factorization"
	print "Super LU (requires ExternalSolversApplication)"
	print "Parallel MKL Pardiso (requires MKLSolversApplication)"
	print "*****************************************************************"
	err
    #else:
	#except LogicError:
    
    if(scaling == False):
        return linear_solver
    else:
        return ScalingSolver(linear_solver,True)

