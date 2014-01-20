from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import benchmarking


def Run():
        print("Running slip condition test...")
        return benchmarking.RunBenchmark("slip_test.py", "slip_test_ref.txt")
