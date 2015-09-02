domain_size = 2


# Fluid solver configuration
class FluidSolverConfiguration:
    solver_type =  "trilinos_vms_monolithic_solver"
    domain_size = 2
    TurbulenceModel = "None"

    # Monolithic solver
    #class linear_solver_config:
    #    solver_type = "Super LU"
    #    scaling = False

    class linear_solver_config:
        solver_type = "AMGCL"
        tolerance = 1E-5
        max_iteration = 5000
        preconditioner = "Diagonal"
        scaling = False
        krylov_type = ""
        smoother_type = ""
    
    #convergence criteria settings
    velocity_relative_tolerance = 1E-3
    velocity_absolute_tolerance = 1E-5
    pressure_relative_tolerance = 1E-3
    pressure_absolute_tolerance = 1E-5
    divergence_cleareance_step = 0
    
    #other solver settings
    oss_switch = 0
    compute_reactions = True
    time_order = 2
    predictor_corrector = False
    dynamic_tau = 0.01
    max_iteration = 10
    laplacian_form = 2

#general problem settings
AutomaticDeltaTime = "Fixed"
Dt = 0.001
Start_time = 0.0
max_time = 1e50
nsteps = 0


groups_dictionary = {
        "fluid" : 5,
                   }
#output settings
output_time = 0.1
output_step = 100
VolumeOutput = True
nodal_results=["NORMAL_SENSITIVITY","SHAPE_SENSITIVITY","ADJOINT_VELOCITY","ADJOINT_PRESSURE","VELOCITY","PRESSURE","REACTION"]
gauss_points_results=[]
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDWriteParticlesFlag = False
GiDMultiFileFlag = "Single"

problem_name="kratosFluid"
problem_path=""
kratos_path="D:\Kratos"
