from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
    Msg = ""
    Text = "========== Adjoint Fluid Application ==========\n"

    #
    # SquareWithUpstreamRectangle

    Text += "SquareWithUpstreamRectangle test: "
    os.chdir("SquareWithUpstreamRectangle/kratos.gid")

    print("running the SquareWithUpstreamRectangle benchmark test ...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "finite_difference_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("SquareWithUpstreamRectangle test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("SquareWithUpstreamRectangle test FAILED")

    os.chdir("../..")

    #
    # Circle

    Text += "Circle test: "
    os.chdir("Circle/kratos.gid")

    print("running the Circle benchmark test ...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "finite_difference_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("Circle test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Circle test FAILED")

    os.chdir("../..")

    #
    # Add other examples here

    #
    print("resume of all of the examples for the AdjointFluidApplication :")
    print(Text)
    return Text

if __name__ == '__main__':
    Run()
