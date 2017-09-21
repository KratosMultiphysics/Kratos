domain_size = 3



class SolverSettings4:
    solver_type = "nonlinear_convection_diffusion_solver"
    domain_size= 3	
    time_order = 1
    predictor_corrector = False
    ReformDofAtEachIteration = False
    echo_level=0
    max_iter = 15;
    toll = 1e-3;

    unknown_variable = "YO"
    density_variable= "DENSITY"
    volume_source_variable= "a"
    diffusion_variable= "rhoD"
    surface_source_variable= "NODAL_PAUX"
    mesh_velocity_variable= "MESH_VELOCITY"
    velocity_variable= "VELOCITY"
    specific_heat_variable= "SPECIFIC_HEAT"
    projection_variable= "TEMP_CONV_PROJ"

    class linear_solver_config: 
        solver_type = "Skyline LU factorization"
        scaling = False           
    

#xsettings2 = ConvectionDiffusionSettings()
#xsettings2.SetDensityVariable(DENSITY)
#xsettings2.SetVolumeSourceVariable(HEAT_FLUX)
#xsettings2.SetSurfaceSourceVariable(NODAL_PAUX)
#xsettings2.SetDiffusionVariable(rhoD)
#xsettings2.SetUnknownVariable(MIXTURE_FRACTION)
#xsettings2.SetMeshVelocityVariable(MESH_VELOCITY)
#xsettings2.SetProjectionVariable(TEMP_CONV_PROJ);
#xsettings2.SetProjectionVariable(THAWTWO);
#xsettings2.SetVelocityVariable(VELOCITY)
#xsettings2.SetSpecificHeatVariable(SPECIFIC_HEAT);
##importing the solvers needed


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



