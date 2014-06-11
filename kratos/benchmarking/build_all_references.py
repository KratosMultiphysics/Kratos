from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
 # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys
import benchmarking

sys.path.append(".")

CurrentDir = os.curdir
KratosRoot = CurrentDir + "/../"
ApplicationsRoot = KratosRoot + "applications/"

os.chdir(ApplicationsRoot)

# SolidMechanicsApplication
os.chdir(ApplicationsRoot + "SolidMechanicsApplication/test_examples")

# PFEM
os.chdir(ApplicationsRoot + "PFEMapplication/test_exemples")
# os.chdir("PFEMapplication/test_examples")

import pfem_build_references

# Fluid Dynamics Application
os.chdir(ApplicationsRoot + "FluidDynamicsApplication/test_examples")

import fluid_dynamics_build_all_references

# Add other directories here
