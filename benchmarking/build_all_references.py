from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
 # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys
import benchmarking

sys.path.append(".")

CurrentDir = os.getcwd()

# SolidMechanicsApplication
os.chdir('../applications')
os.chdir('./SolidMechanicsApplication/test_examples/')


# PFEM
os.chdir('../')
os.chdir('../')
os.chdir('./PFEMapplication/test_exemples')

import pfem_build_references

# Fluid Dynamics Application
os.chdir('../')
os.chdir('../')
os.chdir('./FluidDynamicsApplication/test_examples')

import fluid_dynamics_build_all_references

# Add other directories here
