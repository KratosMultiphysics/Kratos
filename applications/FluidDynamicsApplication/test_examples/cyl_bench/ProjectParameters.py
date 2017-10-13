from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
domain_size = 2

SolverType = "FractionalStep"

# Velocity solver
Velocity_Linear_Solver = "Skyline LU factorization"

# Pressure solver
Pressure_Linear_Solver = "Skyline LU factorization"
TurbulenceModel = "Smagorinsky-Lilly"
SmagorinskyConstant = 0.1

velocity_relative_tolerance = 1E-5
velocity_absolute_tolerance = 1E-6
pressure_relative_tolerance = 1E-5
pressure_absolute_tolerance = 1E-6

time_order = 2
predictor_corrector = False
max_vel_its = 4
max_press_its = 3
laplacian_form = 1

AutomaticDeltaTime = "Fixed"
divergence_cleareance_step = 20
Dt = 0.001
Start_time = 0.0
max_time = 0.1
nsteps = 100

use_dt_in_stabilization = 0.01
use_orthogonal_subscales = 0
Calculate_reactions = True
Use_slip_conditions = True

groups_dictionary = {
    "cyl": 2,
    "vol": 1,
}

output_time = 0.1
output_step = 1  # 100
VolumeOutput = True

nodal_results = ["VELOCITY", "PRESSURE", "REACTION"]
gauss_points_results = []
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDMultiFileFlag = "Single"

problem_name = "cyl_benchFluid"
problem_path = "/home/rrossi/esempi/cyl_bench.gid"
kratos_path = "D:\Kratos"
