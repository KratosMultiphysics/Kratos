from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

kratos_benchmarking_path = '../../../../benchmarking'
sys.path.append(kratos_benchmarking_path)

import benchmarking


def Run():
    print("Running VMS test: parabolic flow in a trapezoidal domain solved using OSS")
    return benchmarking.RunBenchmark("test.py", "trapezoid_exact.txt")

if __name__ == '__main__':
    if Run():
        print("Test successful")
    else:
        print("Test failed")
