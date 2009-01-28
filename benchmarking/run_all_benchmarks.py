import os
import sys
import benchmarking

sys.path.append(".")

Text = "Status of Kratos examples:\n\n"
os.chdir("../applications")

# PFEM
os.chdir("PFEMapplication/test_examples")

import pfem_benchmarks
Text += pfem_benchmarks.Run()

os.chdir("../..")

# Add other directories here

print Text

benchmarking.NotifyViaEmail("Status of Kratos examples", Text, ["pooyan@cimne.upc.edu", "rrossi@cimne.upc.edu", "mossaiby@yahoo.com"])
