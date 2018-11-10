domain_size = 3


# Fluid solver configuration
class FluidSolverConfiguration:
    solver_type =  "vms_monolithic_solver"
    domain_size = 3
    TurbulenceModel = "None"

    # Monolithic solver
    class linear_solver_config:
        solver_type = "AMGCL"
        tolerance = 1E-6
        max_iteration = 5000
        preconditioner = "Diagonal"
        scaling = False

    #convergence criteria settings
    velocity_relative_tolerance = 1E-3
    velocity_absolute_tolerance = 1E-6
    pressure_relative_tolerance = 1E-3
    pressure_absolute_tolerance = 1E-6
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
Dt = 0.16
Start_time = 0.0
max_time = 1
nsteps = 100


groups_dictionary = {
        "fluid" : 1,
                   }
#output settings
output_time = 0.1
output_step = 10
VolumeOutput = True
nodal_results=["VELOCITY","PRESSURE","REACTION"]
gauss_points_results=[]
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDWriteParticlesFlag = False
GiDMultiFileFlag = "Multiples"

problem_name="CandelierFluid"
kratos_path="D:\kratos"
