domain_size = 3

SolverType = "monolithic_solver_eulerian"

class FluidSolverConfiguration:
    solver_type = "SurfaceTension_Monolithic_Solver_3D"
    domain_size = 3
    TurbulenceModel = "None"
    # Monolithic solver
    class linear_solver_config:
        solver_type = "Super LU"
        scaling = False
    
    #convergence criteria settings
    velocity_relative_tolerance = 1E-6
    velocity_absolute_tolerance = 1E-6
    pressure_relative_tolerance = 1E-6
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
Monolithic_Linear_Solver = "BiConjugate gradient stabilized" # "Conjugate gradient" | "MixedUP" | "BiConjugate gradient stabilized"
Monolithic_Iterative_Tolerance = 1E-6
Monolithic_Solver_Max_Iteration = 5000
Monolithic_Preconditioner_type = "ILU0"#"Diagonal"

AutomaticDeltaTime = "Fixed"
divergence_cleareance_step = 10
Dt = 0.00001
Start_time = 0.0
max_time = 10000000000000000000000000.0
nsteps = 10000000000000000000000000

use_dt_in_stabilization = 0.10
use_orthogonal_subscales = 0
Calculate_reactions = True

groups_dictionary = {
        "Fluid" : 1,
                   }

output_time = 0.00001
output_step = 10
VolumeOutput = True

nodal_results=["VELOCITY","PRESSURE","VISCOUS_STRESSX","VISCOUS_STRESSY"]
#nodal_results=["VELOCITY","PRESSURE", "DENSITY", "VISCOSITY",
	       #"IS_INTERFACE", "PHASE_FRACTION_GRADIENT","SOLID_FRACTION_GRADIENT","IS_WATER"]	       
gauss_points_results=[]
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
# Add the following line if using vms_monolithic_solver for lagrangian_model_part
GiDWriteParticlesFlag = False
GiDMultiFileFlag = "Multiples"

problem_name="channel20mm2"
problem_path="/home/alex/Examples_kratos/ALEX/channel30mm6_vms.gid"
kratos_path="home/alex/kratos"
