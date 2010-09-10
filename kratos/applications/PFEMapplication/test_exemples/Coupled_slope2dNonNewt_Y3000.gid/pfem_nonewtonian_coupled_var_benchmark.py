domain_size = 2

Density = 1300.0
Viscosity = 0.0001
Gravity_X = 0.0
Gravity_Y = -9.806
Gravity_Z = 0.0
Porosity = 0.5
Diameter = 0.01
YieldStress = 3000.0
max_x = 2.5
max_y = 3.0
max_z = 1.0
min_x = -1.0
min_y = -1.0
min_z = -1.0
use_oss = 0
dynamic_tau = 1
nsteps = 1000
output_Dt = 0.5
min_dt = 0.01
max_dt = 0.1
safety_factor = 0.1
SolverType = "monolithic_solver_lagrangian"
laplacian_form = 3
fluid_file = r'Coupled_slope2dNonNewt_Y3000_fixed'
# Declare Python Variables

problem_name="Coupled_slope2dNonNewt_Y3000"
problem_path="."
kratos_path="../../../.."
