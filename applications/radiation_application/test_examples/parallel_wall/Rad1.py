domain_size = 2


class SolverSettings2:
    solver_type = "nonlinear_convection_diffusionr_solver"
    domain_size= 3	
    time_order = 1
    predictor_corrector = False
    ReformDofAtEachIteration = False
    echo_level=0
    max_iter = 15;
    toll = 1e-3;

    unknown_variable = "RADIATIVE_INTENSITY_1"
    density_variable= "DENSITY"
    mesh_velocity_variable= "MESH_VELOCITY_1"

    class linear_solver_config: 
        solver_type = "Skyline LU factorization"
        scaling = False           


