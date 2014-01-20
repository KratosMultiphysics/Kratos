from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
kratos_benchmarking_path = '../../../../benchmarking'  # kratos_root/benchmarking

import sys
sys.path.append(kratos_benchmarking_path)

import benchmarking


def Run():
    print("Running remesh.py...")
    return benchmarking.RunBenchmark("remesh.py", "adaptive_mesher2d_ref.txt")
