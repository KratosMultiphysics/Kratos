from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
    Msg = ""
    Text = "== Solid Mechanics ==========\n"

    #
    # VMS2D element test

    Text += "scordelis low roof test: "
    os.chdir("scordelis.gid")
    sys.path.append(os.getcwd())

    print("running the scordelis low roof benchmark test...")
    Msg = benchmarking.RunBenchmark("run_test.py", "min_displacements.txt")

    if (Msg):
        Text += "OK\n"
        print("scordelis low roof test succesful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("scordelis low roof test FAILED")

    os.chdir("..")

    #
    print("resume of all of the examples for the Solid Mechanics application :")
    print(Text)
    return Text

if __name__ == '__main__':
    Run()
