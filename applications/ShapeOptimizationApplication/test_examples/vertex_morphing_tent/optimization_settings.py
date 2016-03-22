class KratosShapeSettings:
    design_history_directory = "./history"
    design_history_file = design_history_directory + "/design_history.csv"
    objectives = {"1":{"grad": "provided"}}
    constraints = {}
    design_control = "vertex_morphing"
    design_output_mode = "relative"
    design_surface_name = "vertex_morphing_tent"
    domain_size = 3
    filter_size = 2.0
    optimization_algorithm = "steepest_descent"
    max_opt_iterations = 200
    relative_tolerance_objective = 1e-3
    normalize_sensitivities = True
    step_size_0 = 0.01
    nodal_results=["OBJECTIVE_SENSITIVITY",
                   "SHAPE_UPDATE",
                   "SHAPE_CHANGE_ABSOLUTE"]
    VolumeOutput = True
    GiDPostMode = "Binary"
    GiDWriteMeshFlag = True
    GiDWriteConditionsFlag = True
    GiDMultiFileFlag = "Single"
