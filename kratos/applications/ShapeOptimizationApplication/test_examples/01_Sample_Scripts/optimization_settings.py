# ===============================================================================================================================
#Optimizer Settings 
# ===============================================================================================================================
class KratosShapeSettings:

    # Workspace settings ---------------------------------------------------------------------------------------------

    design_history_directory = "path/to/desired/Design_history_directory_name"
    design_history_file = design_history_directory + "/design_history.csv"

    # Response functions ---------------------------------------------------------------------------------------------

    # Define container of objective functions
    # Format: objectives = { "unique_func_id": {"grad": "provided"},
    #                        "unique_func_id": {"grad": "provided"},
    #               ... }
    objectives = { "some_objective_name": {"grad": "provided"} }
    # Define container of constraint functions
    # Format: constraints = { "unique_func_id": {"type": "eq"/"ineq", "grad": "provided"},
    #                         "unique_func_id": {"type": "eq"/"ineq", "grad": "provided"},
    #               ... }    
    constraints = { "some_constraint_name": {"type": "eq", "grad": "provided"} }

    # ... general response container possible for tracking of further functions
    
    # Design variables -----------------------------------------------------------------------------------------------

    design_control = "vertex_morphing" 
    # options: "vertex_morphing"

    design_output_mode = "relative"
    # options: "relative"   - X is defined relative to previous design
    #          "total"      - X is defined relative to initial design
    #          "absolute"   - X is defined as absolute values (coordinates)

    # Case: shape_control = "vertex_morphing"
    design_surface_name = "path/to/your/mpda/some_mesh_name"
    domain_size = 3
    filter_size = 0.3

    # Optimization algorithm -----------------------------------------------------------------------------------------

    optimization_algorithm = "steepest_descent" 
    # options: "steepest_descent",
    #          "augmented_lagrange",
    #          "penalized_projection",
    
    # General convergence criterions
    max_opt_iterations = 1000
    
    # Case: "steepest descent"
    relative_tolerance_objective = 1e-1 # [%]
    
    # Case: optimization_algorithm = "augmented_lagrange"
    max_sub_opt_iterations = 100
    relative_tolerance_sub_opt = 1e-1 # [%]
    penalty_fac_0 = 4
    gamma = 4
    penalty_fac_max = 2000
    lambda_0 = 0.0

    # Case: optimization_algorithm = "penalized_projection"
    constraint_scaling = 1000
      
    # Determination of step size (line-search) -----------------------------------------------------------------------

    # Only constant step-size is implemented yet
    step_size_0 = 0.001

    # For GID output -------------------------------------------------------------------------------------------------

    nodal_results=["NORMALIZED_SURFACE_NORMAL",
                   "OBJECTIVE_SENSITIVITY",
                   "CONSTRAINT_SENSITIVITY",
                   "SHAPE_UPDATE",
                   "SHAPE_CHANGE_ABSOLUTE"]
    VolumeOutput = True
    GiDPostMode = "Binary"
    GiDWriteMeshFlag = True
    GiDWriteConditionsFlag = True
    GiDMultiFileFlag = "Single"

    # ----------------------------------------------------------------------------------------------------------------

# ===============================================================================================================================
# Specific parameters of the analyzer
# ===============================================================================================================================

#...

# ===============================================================================================================================