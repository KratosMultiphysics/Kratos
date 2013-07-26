domain_size = 2

SolverType = "FractionalStep"

# Velocity solver
Velocity_Linear_Solver = "Skyline LU factorization"

# Pressure solver
Pressure_Linear_Solver = "Skyline LU factorization"
TurbulenceModel = "None"

velocity_relative_tolerance = 1E-3
velocity_absolute_tolerance = 1E-6
pressure_relative_tolerance = 1E-3
pressure_absolute_tolerance = 1E-6

time_order = 2
predictor_corrector = False
max_vel_its = 4
max_press_its = 3
laplacian_form = 1

Re = 1000

AutomaticDeltaTime = "Fixed"
divergence_cleareance_step = 0
Dt = 1
Start_time = 0.0
max_time = 4
nsteps = 3

use_dt_in_stabilization = 0.01
use_orthogonal_subscales = 0
Calculate_reactions = True

groups_dictionary = {
        "fluid" : 1,
                   }

output_time = 0.1
output_step = 1
VolumeOutput = True

nodal_results=["VELOCITY","PRESSURE","REACTION","DISTANCE"]
gauss_points_results=[]
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDMultiFileFlag = "Single"

problem_name="cavity2d_benchmarkFluid"
problem_path="/home/daniel/kratos/applications/trilinos_application/test_examples/cavity2d_benchmark.gid"
kratos_path="D:\Kratos"
