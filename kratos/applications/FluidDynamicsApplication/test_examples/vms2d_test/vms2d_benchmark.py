from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import benchmarking


def Run():
        print("Running VMS2D test...")
        return benchmarking.RunBenchmark("script_elemtest.py", "vms2d_test_ref.txt")
