import os
import sys

# VMS2D element test

os.chdir("vms2d_test")
sys.path.append(os.getcwd())

import vms2d_build_reference

os.chdir("..")

# slip condition + wall law
os.chdir("slip_test")
sys.path.append(os.getcwd())

import slip_test_build_reference

os.chdir("..")

# Add other examples here

# Note: some test examples are not included in this file, as the results are compared to an exact solution
