##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3

SolverType = "static_poisson_solver"

# Declare Python Variables
	
problem_name="*tcl(file tail [GiD_Info project modelname])"
problem_path="."
#take care: this is for kratos\applications\kMagnetostatic\test_examples\Validation\problem.gid
kratos_path="../../../../../"


