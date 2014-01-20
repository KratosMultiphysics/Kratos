from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

# two_balls_no_damp element test

os.chdir("two_balls_no_damp.gid")
sys.path.append(os.getcwd())

import two_balls_no_damp_build_reference

os.chdir("..")

# two_balls_normal_damp element test

os.chdir("two_balls_normal_damp.gid")
sys.path.append(os.getcwd())

import two_balls_normal_damp_build_reference

os.chdir("..")

# rotating_ball_no_tangent_damp element test

os.chdir("rotating_ball_no_tangent_damp.gid")
sys.path.append(os.getcwd())

import rotating_ball_no_tangent_damp_build_reference

os.chdir("..")

# rotating_ball_rolling_friction element test

os.chdir("rotating_ball_rolling_friction.gid")
sys.path.append(os.getcwd())

import rotating_ball_rolling_friction_build_reference

os.chdir("..")
# Add other examples here

# Note: some test examples are not included in this file, as the results are compared to an exact solution
