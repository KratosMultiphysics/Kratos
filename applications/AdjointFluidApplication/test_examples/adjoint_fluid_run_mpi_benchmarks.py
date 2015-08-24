from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
    Msg = ""
    Text = "========== MPI Adjoint Fluid Application ==========\n"

    #
    # NACA0012

    Text += "NACA0012 test: "
    os.chdir("NACA0012/kratos.gid")

    print("running the NACA0012 benchmark test ...")
    successful,Msg = benchmarking.MPIParallelRunBenchmark("run_test.py", "finite_difference_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("NACA0012 test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("NACA0012 test FAILED")

    os.chdir("../..")

    #
    # Add other examples here

    #
    print("resume of all of the examples for the MPI AdjointFluidApplication :")
    print(Text)
    return Text

if __name__ == '__main__':
    Run()
