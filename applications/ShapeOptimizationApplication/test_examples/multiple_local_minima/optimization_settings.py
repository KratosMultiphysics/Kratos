design_history_directory = "Design_history"
design_history_file = "design_history.csv"
objectives = {"1":{"grad": "provided"}}
constraints = {}
design_control = "vertex_morphing"
design_output_mode = "relative"
design_surface_name = "multiple_local_minima"
domain_size = 3
filter_function = "linear"
use_mesh_preserving_filter_matrix = False
filter_size = 5.0 # 7.0
perform_edge_damping = False
damped_edges = []
optimization_algorithm = "steepest_descent"
max_opt_iterations = 200
relative_tolerance_objective = 1e-3
normalize_search_direction = True
step_size = 0.01
nodal_results=["OBJECTIVE_SENSITIVITY",
               "SHAPE_UPDATE",
               "SHAPE_CHANGE_ABSOLUTE"]
VolumeOutput = True
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDMultiFileFlag = "Single"
