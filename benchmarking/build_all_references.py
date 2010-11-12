import os
import sys
import benchmarking

sys.path.append(".")

CurrentDir = os.curdir
KratosRoot = CurrentDir + "/../"
ApplicationsRoot = KratosRoot + "applications/"

os.chdir(ApplicationsRoot)

# PFEM
os.chdir(ApplicationsRoot + "PFEMapplication/test_exemples")
#os.chdir("PFEMapplication/test_examples")

import pfem_build_references

# Fluid Dynamics Application
os.chdir(ApplicationsRoot + "FluidDynamicsApplication/test_examples")

import fluid_dynamics_build_all_references

# Add other directories here
