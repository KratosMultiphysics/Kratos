from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import os
#Declaring Problem name and path
#Declaring Problem name and path
current_directory=os.getcwd()
problem_path = current_directory
os.chdir("..\\..\\..\\..\\..")
kratos_path = os.getcwd()
os.chdir(current_directory)
#os.chdir(current_directory)

# including kratos path
kratos_benchmarking_path=os.path.join(kratos_path,'kratos\\benchmarking')
sys.path.append(kratos_benchmarking_path)

import benchmarking


def Run():
        print("Running  SUPGConvDiffPhaseChange3D - CUBE test.")
        return benchmarking.RunBenchmark("cube_benchmark.py", "cube_benchmark_test_ref.txt")
