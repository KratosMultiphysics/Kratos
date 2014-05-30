from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
 # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys
import benchmarking

sys.path.append(".")

CurrentDir = os.getcwd()
KratosRoot = CurrentDir + "/../"
ApplicationsRoot = KratosRoot + "applications/"

Text = "Status of Kratos examples:\n\n"
os.chdir(ApplicationsRoot)

# SolidMechanics application
os.chdir(ApplicationsRoot + "SolidMechanicsApplication/test_examples")
import solid_mechanics_run_all_benchmarks
Text += solid_mechanics_run_all_benchmarks.Run()

# Fluid Dynamics application
os.chdir(ApplicationsRoot + "FluidDynamicsApplication/test_examples")
import fluid_dynamics_run_all_benchmarks
Text += fluid_dynamics_run_all_benchmarks.Run()

# DEM application
os.chdir(ApplicationsRoot + "DEM_application/test_examples")
import DEM_benchmarks
Text += DEM_benchmarks.Run()

# Structural application
os.chdir(ApplicationsRoot + "structural_application/test_examples")
import structural_benchmarks
Text += structural_benchmarks.Run()

# Incompressible fluid application
os.chdir(ApplicationsRoot + "incompressible_fluid_application/test_examples")
import incompressible_benchmarks
Text += incompressible_benchmarks.Run()

# Convection-Diffusion application
os.chdir(ApplicationsRoot + "convection_diffusion_application/test_examples")
import convection_diffusion_benchmarks
Text += convection_diffusion_benchmarks.Run()

# FSI application
os.chdir(ApplicationsRoot + "FSIapplication/test_examples")
import FSI_benchmarks
Text += FSI_benchmarks.Run()

# PFEM application
# os.chdir(ApplicationsRoot + "PFEMapplication/test_exemples")
# import pfem_benchmarks
# Text += pfem_benchmarks.Run()

# Thermo Mechanical application
os.chdir(ApplicationsRoot + "ThermoMechanicalApplication/test_examples")
import ThermoMechanical_benchmarks
Text += ThermoMechanical_benchmarks.Run()

# Meshing application
os.chdir(ApplicationsRoot + "MeshingApplication/test_exemples")
import meshing_benchmarks
Text += meshing_benchmarks.Run()

# Trilinos application
os.chdir(ApplicationsRoot + "trilinos_application/test_examples")
import trilinos_benchmarks
Text += trilinos_benchmarks.Run()

# Add other directories here

# print Text
os.chdir(CurrentDir)
outputfile = open("benchmarking_output.txt", "w")
outputfile.write(Text)


print(Text)

# benchmarking.NotifyViaEmail("Status of Kratos examples", Text, ["pooyan@cimne.upc.edu", "rrossi@cimne.upc.edu", "mossaiby@yahoo.com"])
