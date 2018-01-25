from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
domain_size = 3


class SolverSettings1:
    solver_type = "convection_diffusion_solver"
    domain_size = 2
    time_order = 1
    predictor_corrector = False
    ReformDofAtEachIteration = True
    echo_level = 0

    # set variables to be used
    unknown_variable = "TEMPERATURE"
    density_variable = "DENSITY"
    projection_variable = "TEMP_CONV_PROJ"
    volume_source_variable = "HEAT_FLUX"
    diffusion_variable = "CONDUCTIVITY"
    surface_source_variable = "FACE_HEAT_FLUX"
    mesh_velocity_variable = "MESH_VELOCITY"
    velocity_variable = "VELOCITY"
    specific_heat_variable = "SPECIFIC_HEAT"

    class linear_solver_config:
        solver_type = "Conjugate gradient"
        max_iteration=5000
        tolerance=0.0001
        scaling = False


class SolverSettings2:
    solver_type = "nonlinear_convection_diffusion_solver"
    domain_size = 2
    time_order = 1
    predictor_corrector = False
    ReformDofAtEachIteration = True
    echo_level = 0
    max_iter = 15
    toll = 1e-3

    # set variables to be used
    unknown_variable = "TEMPERATURE"
    density_variable = "DENSITY"
    projection_variable = "TEMP_CONV_PROJ"
    volume_source_variable = "HEAT_FLUX"
    diffusion_variable = "CONDUCTIVITY"
    surface_source_variable = "FACE_HEAT_FLUX"
    mesh_velocity_variable = "MESH_VELOCITY"
    velocity_variable = "VELOCITY"
    specific_heat_variable = "SPECIFIC_HEAT"

    class linear_solver_config:
        solver_type = "Skyline LU factorization"
        scaling = False


Dt = 0.1
Start_time = 0.0
max_time = 60
max_time = 40
nsteps = 300
output_step = 1
output_time = 0.0

nodal_results = ["TEMPERATURE"]
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDMultiFileFlag = "Single"

#problem_name = "square"
#problem_path = "/home/julio.marti/new_kratos/applications/convection_diffusion_application/test_examples/square.gid"
kratos_path = "../../../.."
