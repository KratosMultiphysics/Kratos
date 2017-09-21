domain_size = 3



class SolverSettings3:
    solver_type = "nonlinear_convection_diffusion_solver"
    domain_size= 3	
    time_order = 1
    predictor_corrector = False
    ReformDofAtEachIteration = False
    echo_level=0
    max_iter = 15;
    toll = 1e-3;

    unknown_variable = "YF"
    density_variable= "DENSITY"
    volume_source_variable= "a"
    diffusion_variable= "M"
    surface_source_variable= "NODAL_PAUX"
    mesh_velocity_variable= "MESH_VELOCITY"
    velocity_variable= "VELOCITY"
    specific_heat_variable= "SPECIFIC_HEAT"
    projection_variable= "TEMP_CONV_PROJ"

    class linear_solver_config: 
        solver_type = "Skyline LU factorization"
        scaling = False           
    



Dt = 0.1
Start_time = 0.0
max_time = 60
max_time = 40
nsteps = 300
output_step=1
output_time = 0.0

#nodal_results=["TEMPERATURE"]
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDMultiFileFlag = "Single"

#problem_name="square"
#problem_path="/home/julio.marti/new_kratos/applications/convection_diffusion_application/test_examples/square.gid"
#kratos_path="../../../.."



