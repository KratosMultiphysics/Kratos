from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import benchmarking


def Run():
    print("Running dam2d.py...")
    return benchmarking.RunBenchmark("dam2d.py", "dam2d_ref.txt")
