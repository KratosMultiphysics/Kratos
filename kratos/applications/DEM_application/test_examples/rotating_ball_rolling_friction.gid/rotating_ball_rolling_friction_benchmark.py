from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import benchmarking


def Run():
    print("Running rotating_ball_rolling_friction test...")
    return benchmarking.RunBenchmark("rotating_ball_rolling_friction.py", "rotating_ball_rolling_friction_ref.txt")
