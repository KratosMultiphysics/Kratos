from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
kratos_benchmarking_path = '../../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking

print("running reference data for edgebased_PureConvection.py...")
benchmarking.RunBenchmark("edgebased_PureConvection.py", "test_pureconvectionsolver_benchmarking_ref.txt")
