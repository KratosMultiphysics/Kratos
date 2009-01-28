import os
import sys
import benchmarking

sys.path.append(".")

os.chdir("../applications")

# PFEM
os.chdir("PFEMapplication/test_examples")

import pfem_build_references

os.chdir("../..")

# Add other directories here
