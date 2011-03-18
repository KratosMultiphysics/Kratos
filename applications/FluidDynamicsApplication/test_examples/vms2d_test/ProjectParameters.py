domain_size = 2

SolverType = "pressure_splitting"
#SolverType = "monolithic_solver_eulerian"
#Monolithic_Linear_Solver = "Skyline LU factorization"
Monolithic_Linear_Solver = "BiConjugate gradient stabilized"
Monolithic_Iterative_Tolerance = 1E-6
Monolithic_Solver_Max_Iteration = 5000
Monolithic_Preconditioner_type = "Diagonal"
# Velocity solver
#Velocity_Linear_Solver = "Super LU"
Velocity_Linear_Solver = "BiConjugate gradient stabilized"
Velocity_Iterative_Tolerance = 1E-6
Velocity_Solver_Max_Iteration = 5000
Velocity_Preconditioner_type = "Diagonal"

# Pressure solver
#Pressure_Linear_Solver = "Super LU"
Pressure_Linear_Solver = "BiConjugate gradient stabilized"
Pressure_Iterative_Tolerance = 1E-6
Pressure_Solver_Max_Iteration = 5000
Pressure_Preconditioner_type = "Diagonal"

velocity_relative_tolerance = 1E-3
velocity_absolute_tolerance = 1E-6
pressure_relative_tolerance = 1E-3
pressure_absolute_tolerance = 1E-6

time_order = 2
predictor_corrector = False
max_iterations = 10
laplacian_form = 2

Dt = 0.05
Start_time = 0.0
##max_time = 60
max_time = 40
nsteps = 1200

use_dt_in_stabilization = 0.01
use_orthogonal_subscales = 0
Calculate_reactions = True

output_time = 0.0
output_step = 0

nodal_results=["VELOCITY","PRESSURE","REACTION","ADVPROJ","DIVPROJ"]
gauss_points_results=["VORTICITY"]
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDMultiFileFlag = "Single"

problem_name="vms2d_testFluid"
problem_path="/home/jcotela/Problems/myelem_test.gid"
kratos_path="../../../.."
