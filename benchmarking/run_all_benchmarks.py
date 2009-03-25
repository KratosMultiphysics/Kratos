import os
import sys
import benchmarking

sys.path.append(".")

CurrentDir = os.curdir
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


# Add other directories here

print Text

#benchmarking.NotifyViaEmail("Status of Kratos examples", Text, ["pooyan@cimne.upc.edu", "rrossi@cimne.upc.edu", "mossaiby@yahoo.com"])
