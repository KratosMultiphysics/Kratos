def import_solver( SolverSettings ):
	"""this function imports a solver named "solver_type" from SolverSettings
	solver_type is expected to be the FILENAME of the solver to be imported"""
	obj =  __import__(SolverSettings.solver_type)
	return obj
