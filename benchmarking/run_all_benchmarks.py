import os
import sys
import benchmarking

sys.path.append(".")

CurrentDir = os.getcwd()
KratosRoot = CurrentDir + "/../"
ApplicationsRoot = KratosRoot + "applications/"

Text = "Status of Kratos examples:\n\n"
os.chdir(ApplicationsRoot)

# Structural application                                                                                                                                                                           
os.chdir(ApplicationsRoot + "structural_application/test_examples")
import structural_benchmarks
Text += structural_benchmarks.Run()


# Incompressible fluid application
os.chdir(ApplicationsRoot + "incompressible_fluid_application/test_examples")
import incompressible_benchmarks
Text += incompressible_benchmarks.Run()


# Meshing application
os.chdir(ApplicationsRoot + "MeshingApplication/test_exemples")
import meshing_benchmarks
Text += meshing_benchmarks.Run()

# Convection-Diffusion application
os.chdir(ApplicationsRoot + "convection_diffusion_application/test_examples")
import convection_diffusion_benchmarks
Text += convection_diffusion_benchmarks.Run()

# FSI application
os.chdir(ApplicationsRoot + "FSIapplication/test_examples")
import FSI_benchmarks
Text += FSI_benchmarks.Run()

### PFEM application
##os.chdir(ApplicationsRoot + "PFEMapplication/test_exemples")
##import pfem_benchmarks
##Text += pfem_benchmarks.Run()

# Fluid Dynamics application
os.chdir(ApplicationsRoot + "FluidDynamicsApplication/test_examples")
import fluid_dynamics_run_all_benchmarks
Text += fluid_dynamics_run_all_benchmarks.Run()

# Add other directories here

#print Text
os.chdir(CurrentDir)
outputfile = open("benchmarking_output.txt","w")
outputfile.write(Text)


print Text

#benchmarking.NotifyViaEmail("Status of Kratos examples", Text, ["pooyan@cimne.upc.edu", "rrossi@cimne.upc.edu", "mossaiby@yahoo.com"])
