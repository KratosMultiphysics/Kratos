from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
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

# fractional step: 2d cavity
os.chdir("fs_cavity")
sys.path.append(os.getcwd())

import fs_cavity_build_reference

os.chdir("..")

# test of cylinder with slip BC, wall law and outlet close to the cylinder
os.chdir("cyl_bench")
sys.path.append(os.getcwd())

import cyl_bench_build_reference

os.chdir("..")

# Add other examples here

# Note: some test examples are not included in this file, as the results are compared to an exact solution
