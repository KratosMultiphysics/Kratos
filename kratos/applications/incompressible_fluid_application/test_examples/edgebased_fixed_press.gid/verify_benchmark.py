from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
kratos_benchmarking_path = '../../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking

print("running the benchmark for edgebased_fixed_press test...")
successful,Msg = benchmarking.RunBenchmark(
    "run_benchmark.py",
    "benchmark_reference_solution.txt")

if(successful==True):
    print("edgebased_fixed_press test example successful")
else:
    print("edgebased_fixed_press test example FAILED")
