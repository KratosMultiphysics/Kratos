






# "algorithm_settings": {
#         "type": "algorithm_relaxed_gradient_projection",
#         "settings": {
#             "echo_level": 0,
#             "line_search": {
#                 "type": "const_step",
#                 "init_step": 1e-2,
#                 "gradient_scaling": "inf_norm"
#             },
#             "conv_settings": {
#                 "type": "max_iter",
#                 "max_iter": 10
#             },
#             "linear_solver_settings": {
#                 "solver_type": "LinearSolversApplication.dense_col_piv_householder_qr"
#             }
#         },
#         "controls": [
#             "thickness_control"
#         ],
#         "objective": {
#             "response_name": "mass",
#             "type": "minimization",
#             "scaling": 1.0
#         },
#         "constraints": [
#             {
#                 "response_name": "strain_energy",
#                 "type": "<=",
#                 "scaled_ref_value": "initial_value"
#             }
#         ]
#     }