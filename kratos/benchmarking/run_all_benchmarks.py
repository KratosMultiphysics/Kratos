import os
import sys
import benchmarking

sys.path.append(".")

CurrentDir = os.curdir
KratosRoot = CurrentDir + "/../"
ApplicationsRoot = KratosRoot + "applications/"

Text = "Status of Kratos examples:\n\n"
os.chdir(ApplicationsRoot)

# PFEM                                                                                                                                                                           
os.chdir(ApplicationsRoot + "PFEMapplication/test_exemples")

import pfem_benchmarks
Text += pfem_benchmarks.Run()

# Add other directories here

print Text

benchmarking.NotifyViaEmail("Status of Kratos examples", Text, ["pooyan@cimne.upc.edu", "rrossi@cimne.upc.edu", "mossaiby@yahoo.com"])
