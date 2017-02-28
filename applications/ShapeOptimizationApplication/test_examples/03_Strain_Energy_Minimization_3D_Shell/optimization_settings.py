# ================================================================================================================
# Workspace settings 
# ================================================================================================================

# Define directory with a relative or absolute path
design_history_directory = "Design_history"
design_history_file = "design_history.csv"

# ================================================================================================================
# Response functions 
# ================================================================================================================

# Define container of objective functions
# Format: objectives = { "unique_func_id": {"gradient_mode": "analytic"},
#                        "unique_func_id": {"gradient_mode": "semi_analytic", "step_size": 1e-5},
#                        "unique_func_id": {"gradient_mode": "external"},
#                        ... }
objectives = { "strain_energy": {"gradient_mode": "semi_analytic", "step_size": 1e-8} }

# Define container of constraint functions
# Format: constraints = { "unique_func_id": {"type": "eq"/"ineq","gradient_mode": "analytic"},
#                         "unique_func_id": {"type": "eq"/"ineq","gradient_mode": "semi_analytic", "step_size": 1e-5},
#                         "unique_func_id": {"type": "eq"/"ineq","gradient_mode": "external"},
#                         ... }    
constraints = {  }
    
# ================================================================================================================  
# Design variables 
# ================================================================================================================

design_control = "vertex_morphing" 
# options: "vertex_morphing"

design_output_mode = "relative"
# options: "relative"   - X is defined relative to previous design
#          "total"      - X is defined relative to initial design
#          "absolute"   - X is defined as absolute values (coordinates)


# ================================================================================================================
# Case: design_control = "vertex_morphing"
# ================================================================================================================

model_input_filename = "3D_Shell"
design_surface_sub_model_name = "design_surface"
domain_size = 3
# options: 2 or 3 for 2D or 3D optimization patch

filter_function = "linear"
# options: "gaussian"
#          "linear"
filter_size = 3

use_mesh_preserving_filter_matrix = False
# options: True    - surface normal information used in the filter matrix
#        : False   - complete filter matrix is used

perform_damping = True
# options: True    - damping is applied with the settings below
#        : False   - no damping is applied, settings below can be ignored

damping_regions = [ ["support_edges", False, True, True, "linear", 3 ],
                    ["side_edges", False, False, True, "linear", 3 ] ]
# damping_region = [ [sub_model_part_name_1, damp_in_X, damp_in_Y, damp_in_Z, damping_function, damping_radius ],
#                    [sub_model_part_name_2, damp_in_X, damp_in_Y, damp_in_Z, damping_function, damping_radius ],
#                    ... ]
# options for damping function: "cosine"
#                               "linear"

# ================================================================================================================
# Optimization algorithm 
# ================================================================================================================

optimization_algorithm = "steepest_descent" 
# options: "steepest_descent",
#          "augmented_lagrange",
#          "penalized_projection",
    
# General convergence criterions
max_opt_iterations = 300
    
# Case: "steepest descent"
relative_tolerance_objective = 1e-1 # [%]
    
# Case: optimization_algorithm = "augmented_lagrange"
max_sub_opt_iterations = 100
relative_tolerance_sub_opt = 1e-1 # [%]
penalty_fac_0 = 4
gamma = 4
penalty_fac_max = 2000
lambda_0 = 0.0
            
# ================================================================================================================ 
# Determination of step size (line-search) 
# ================================================================================================================

# Only constant step-size is implemented yet
normalize_search_direction = True
step_size = .1 # e.g. 5 for active normalization or 1e7 for inactive normalization

# ================================================================================================================
# For GID output 
# ================================================================================================================

nodal_results=[ "NORMALIZED_SURFACE_NORMAL",
                "OBJECTIVE_SENSITIVITY",
                "MAPPED_OBJECTIVE_SENSITIVITY",
                "DESIGN_UPDATE",
                "DESIGN_CHANGE_ABSOLUTE",
                "SHAPE_UPDATE",
                "SENSITIVITIES_DEACTIVATED",
                "SHAPE_UPDATES_DEACTIVATED",
                "SHAPE_CHANGE_ABSOLUTE"]
VolumeOutput = True
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDMultiFileFlag = "Single"

# ================================================================================================================
