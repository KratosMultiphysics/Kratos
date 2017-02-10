design_history_directory = "Design_history"
design_history_file = "design_history.csv"
objectives = {"1":{"grad": "provided"}}
constraints = {}
design_control = "vertex_morphing"
design_output_mode = "relative"
design_surface_name = "vertex_morphing_tent"
domain_size = 3
filter_function = "gaussian"
filter_size = 2.0
use_mesh_preserving_filter_matrix = False
perform_edge_damping = False
damped_edges = []
optimization_algorithm = "steepest_descent"
max_opt_iterations = 200
relative_tolerance_objective = 1e-3
normalize_search_direction = False
step_size = 0.01
nodal_results=["OBJECTIVE_SENSITIVITY",
               "SHAPE_UPDATE",
               "SHAPE_CHANGE_ABSOLUTE"]
VolumeOutput = True
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDMultiFileFlag = "Single"
