from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

#Declaring Problem name and path
kratos_path = "\\KratosMultiphysics"

# including kratos path
kratos_benchmarking_path=kratos_path+'\\kratos\\benchmarking'
sys.path.append(kratos_path)
sys.path.append(kratos_benchmarking_path)
from KratosMultiphysics import *
import benchmarking

print("Building reference data for SUPGConvDiffPhaseChange3D - CUBE test...")
benchmarking.BuildReferenceData("cube_benchmark.py", "cube_benchmark_test_ref.txt")