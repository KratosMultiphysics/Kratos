domain_size = 2

SolverType = "monolithic_solver_eulerian"
#SolverType2 = "FractionalStep"

class FluidSolverConfiguration:
    solver_type =  "SurfaceTension_monolithic_solver"
    domain_size = 2
    TurbulenceModel = "None"

    # Monolithic solver
    class linear_solver_config:
        solver_type = "Super LU"
        scaling = False
    #convergence criteria settings
    velocity_relative_tolerance = 1E-4
    velocity_absolute_tolerance = 1E-6
    pressure_relative_tolerance = 1E-4
    pressure_absolute_tolerance = 1E-6
    divergence_cleareance_step = 1
    #other solver settings
    oss_switch = 0
    compute_reactions = True
    time_order = 2
    predictor_corrector = False
    dynamic_tau = 0.1
    max_iteration = 20
    laplacian_form = 2
    eulerian_model_part = 0
# Monolithic solver
Monolithic_Linear_Solver ="MixedUP"#"BiConjugate gradient stabilized"#
Monolithic_Iterative_Tolerance = 1E-4
Monolithic_Solver_Max_Iteration = 5000
Monolithic_Preconditioner_type = "ILU0"#"Diagonal"
Velocity_Linear_Solver="BiConjugate gradient stabilized"
Pressure_Linear_Solver="Conjugate gradient"
Velocity_Preconditioner_type="ILU0"
Pressure_Preconditioner_type="ILU0"
Velocity_Iterative_Tolerance=1E-6
Pressure_Iterative_Tolerance=1E-3
Velocity_Solver_Max_Iteration = 5000
Pressure_Solver_Max_Iteration = 1000

TurbulenceModel = "None"

velocity_relative_tolerance = 1E-4
velocity_absolute_tolerance = 1E-6
pressure_relative_tolerance = 1E-4
pressure_absolute_tolerance = 1E-6

time_order = 2
predictor_corrector = False
max_iterations = 10
laplacian_form = 2

AutomaticDeltaTime = "Fixed"
divergence_cleareance_step = 10
Dt = 0.01
Start_time = 0.0
max_time = 1.00
nsteps = 100

use_dt_in_stabilization = 0.10
use_orthogonal_subscales = 0
Calculate_reactions = True

groups_dictionary = {
        "Fluid" : 1,
                   }

output_time = 0.01
output_step = 10
VolumeOutput = True

nodal_results=["VELOCITY","PRESSURE"]
gauss_points_results=[]
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
# Add the following line if using SurfaceTension_monolithic_solver for lagrangian_model_part
GiDWriteParticlesFlag = False
GiDMultiFileFlag = "Multiples"


kratos_path="home/alex/kratos"
