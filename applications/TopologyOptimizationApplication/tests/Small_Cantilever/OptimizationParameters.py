#========================================================================================
# RESPONSE FUNCTIONS
#========================================================================================

# Define container of objective functions
# Format: objectives = { "unique_func_id": {"grad": "provided"},
#                        "unique_func_id": {"grad": "provided"},
#               		... }
objectives = { "strain_energy": {"grad": "provided"} }

# Define container of constraint functions
# Format: constraints = { "unique_func_id": {"type": "eq"/"ineq", "grad": "provided"},
#                         "unique_func_id": {"type": "eq"/"ineq", "grad": "provided"},
#               		... }    
constraints = { "volume_fraction": {"type": "eq", "grad": "provided"} }

#========================================================================================
