from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
domain_size = 3

max_time = 2.00000e+00
safety_factor = 5.00000e-01
max_time_step = 1.00000e-03
output_dt = 1.00000e-01
density = 1.00000e+03
viscosity = 1.00000e-06
body_force_x = 0.00000e+00
body_force_y = 0.00000e+00
body_force_z = -9.81000e+00
wall_law_y = 0.00000e+00
use_mass_correction = True
redistance_frequency = 5.00000e+00
extrapolation_layers = 5.00000e+00
number_of_inital_steps = 1.00000e+01
initial_time_step = 1.00000e-05
reduction_on_failure = 3.00000e-01
stabdt_pressure_factor = 1.00000e+00
stabdt_convection_factor = 1.00000e-02
tau2_factor = 1.00000e+00
edge_detection_angle = 4.50000e+01
assume_constant_pressure = False
compute_porous_resistance_law = 0.00000e+00
pressure_fixed = "ON"
fix_location = "inlet"
Z_coord_free_surface = 5.20000e-01
Z_coord_bottom = 0.00000e+00
X_coord_I_O = 0.00000e+00
free_surface_Z_coord = 1.00000e+00
X1 = 7.50000e-01
Y1 = 2.50000e-01
X2 = 1.00000e+00
Y2 = 1.00000e+00
X3 = 1.00000e+00
Y3 = 1.00000e+00
X4 = 1.00000e+00
Y4 = 1.00000e+00
X5 = 1.00000e+00
Y5 = 1.00000e+00
Xtol = 5.00000e-02
Ytol = 5.00000e-02
# Declare Python Variables

problem_name = 'edgebased_fixed_press'
problem_path = '/home/CIMNE/Desktop/edgebased_fixed_press.gid'
kratos_path = '/home/CIMNE/kratos'
