import os
import sys

# two_balls_no_damp element test

os.chdir("two_balls_normal_damp.gid")
sys.path.append(os.getcwd())

import two_balls_normal_damp_build_reference

os.chdir("..")

# Add other examples here

# Note: some test examples are not included in this file, as the results are compared to an exact solution
