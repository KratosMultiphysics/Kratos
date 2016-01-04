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



# Kratos core tests
try:
    os.chdir(KratosRoot + "kratos/test_examples")
    import kratos_run_all_benchmarks
    Text += kratos_run_all_benchmarks.Run()
except:
    pass


os.chdir(ApplicationsRoot)

# SolidMechanics application
try:
    os.chdir(ApplicationsRoot + "SolidMechanicsApplication/test_examples")
    import solid_mechanics_run_all_benchmarks
    Text += solid_mechanics_run_all_benchmarks.Run()
except:
    pass

# Fluid Dynamics application
try:
    os.chdir(ApplicationsRoot + "FluidDynamicsApplication/test_examples")
    import fluid_dynamics_run_all_benchmarks
    Text += fluid_dynamics_run_all_benchmarks.Run()
except:
    pass    

# DEM application
try:
    os.chdir(ApplicationsRoot + "DEM_application/test_examples")
    import DEM_run_all_benchmarks
    Text += DEM_run_all_benchmarks.Run()
except:
    pass    

# Structural application
try:
    os.chdir(ApplicationsRoot + "structural_application/test_examples")
    import structural_benchmarks
    Text += structural_benchmarks.Run()
except:
    pass    

# Incompressible fluid application
try:
    os.chdir(ApplicationsRoot + "incompressible_fluid_application/test_examples")
    import incompressible_benchmarks
    Text += incompressible_benchmarks.Run()
except:
    pass    

# Convection-Diffusion application
try:
    os.chdir(ApplicationsRoot + "convection_diffusion_application/test_examples")
    import convection_diffusion_benchmarks
    Text += convection_diffusion_benchmarks.Run()
except:
    pass    

# FSI application
try:
    os.chdir(ApplicationsRoot + "FSIapplication/test_examples")
    import FSI_benchmarks
    Text += FSI_benchmarks.Run()
except:
    pass    

# PFEM application
#try:
#   os.chdir(ApplicationsRoot + "PFEMapplication/test_exemples")
#   import pfem_benchmarks
#   Text += pfem_benchmarks.Run()
#except:
#    pass

# Thermo Mechanical application
try:
    os.chdir(ApplicationsRoot + "ThermoMechanicalApplication/test_examples")
    import ThermoMechanical_benchmarks
    Text += ThermoMechanical_benchmarks.Run()
except:
    pass

# Meshing application
try:
    os.chdir(ApplicationsRoot + "MeshingApplication/test_exemples")
    import meshing_benchmarks
    Text += meshing_benchmarks.Run()
except:
    pass    

# Trilinos application
try:
    os.chdir(ApplicationsRoot + "trilinos_application/test_examples")
    import trilinos_benchmarks
    Text += trilinos_benchmarks.Run()
except:
    pass    

# SwimmingDem application
try:
    os.chdir(ApplicationsRoot + "swimming_DEM_application/test_examples")
    import swimming_DEM_run_all_benchmarks
    Text += swimming_DEM_run_all_benchmarks.Run()
except:
    pass

# Add other directories here

# print Text
os.chdir(CurrentDir)
outputfile = open("benchmarking_output.txt", "w")
outputfile.write(Text)


print(Text)

# benchmarking.NotifyViaEmail("Status of Kratos examples", Text, ["pooyan@cimne.upc.edu", "rrossi@cimne.upc.edu", "mossaiby@yahoo.com"])
