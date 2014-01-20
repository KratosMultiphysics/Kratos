from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# -*- coding: utf-8 -*-
domain_size = 2

Dt = 0.1
max_time = 10.0
SolverType = "fractional_step"
laplacian_form = 1
output_step = 10
predictor_corrector = True
max_press_its = 3
use_oss = 0
dynamic_tau = 0.001
# Declare Python Variables

problem_name = "cavity2D"
problem_path = "."
kratos_path = "../../../.."
