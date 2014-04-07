from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
    Msg = ""
    Text = "== ThermoMechanical Application ==========\n"

    #
    # naca3d test for edge based level set solver

    Text += "Activation_test : "
    os.chdir("Activation_test.gid")
    sys.path.append(os.getcwd())

    print("running the benchmark for the Activation_test ...")
    successful,Msg = benchmarking.RunBenchmark("run_benchmark.py", "Activation_reference_solution.txt")

    if(successful==True):
        Text += "OK\n"
        print("Activation_test  example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Activation_test example FAILED")

    os.chdir("..")

    #

    # Add other examples here

    # Add other examples here
    #
    print("resume of all of the examples for the ThermoMechanical Application  :")
    print(Text)
    return Text

Run()
