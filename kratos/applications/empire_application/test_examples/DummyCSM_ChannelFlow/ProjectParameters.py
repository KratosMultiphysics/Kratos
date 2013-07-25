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

AutomaticDeltaTime = "Fixed"
divergence_cleareance_step = 0
Dt = 0.001
Start_time = 0.0
max_time = 1
nsteps = 1000

use_dt_in_stabilization = 0.01
use_orthogonal_subscales = 0
Calculate_reactions = True

groups_dictionary = {
        "fluid" : 1,
                   }

output_time = 0.1
output_step = 100
VolumeOutput = True

nodal_results=["VELOCITY","PRESSURE","REACTION","DISTANCE"]
gauss_points_results=[]
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDMultiFileFlag = "Single"

problem_name="ChannelFlow"
problem_path="/home/jwolf/kratos/applications/empire_application/test_examples/ChannelFlow"
kratos_path="/home/jwolf/kratos"
