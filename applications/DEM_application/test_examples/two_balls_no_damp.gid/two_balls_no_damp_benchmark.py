from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import benchmarking


def Run():
    print("Running two_balls_no_damp test...")
    return benchmarking.RunBenchmark("two_balls_no_damp.py", "two_balls_no_damp_ref.txt")
