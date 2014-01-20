from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
domain_size = 3

Dt = 0.02
max_time = 2.0
SolverType = "fractional_step"
laplacian_form = 1
output_step = 10
predictor_corrector = False
max_press_its = 3
use_oss = 0
dynamic_tau = 0.0  # 1.0
# Declare Python Variables

problem_name = "naca3d"
problem_path = "/home/rrossi/esempi/naca3d.gid"
kratos_path = "/home/rrossi/kratos"
