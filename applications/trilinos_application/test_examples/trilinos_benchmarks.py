from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
    Msg = ""
    Text = "========== Trilinos Aplication ==========\n"

    #

    Text += "Trilinos -- Cantilever_Amesos: "
    os.chdir("cantilever2d.gid")
    sys.path.append(os.getcwd())

    print("Running Trilinos -- Cantilever_Amesos...")
    successful, Msg = benchmarking.MPIParallelRunBenchmark("cantilever_amesos.py", "cantilever.txt")

    if(successful==True):
        Text += "OK\n"
        print("Trilinos Cantilever_Amesos example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("cantilever_amesos.py example FAILED")

    os.chdir("..")

    #
    Text += "Trilinos -- Cantilever_Aztec: "
    os.chdir("cantilever2d.gid")
    sys.path.append(os.getcwd())

    print("Running Trilinos -- Cantilever_Aztec...")
    success, Msg = benchmarking.MPIParallelRunBenchmark("cantilever_aztec.py", "cantilever.txt")

    if(successful==True):
        Text += "OK\n"
        print("Trilinos Cantilever_Aztec example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("cantilever_aztec.py example FAILED")

    os.chdir("..")

    #
    Text += "Trilinos -- Cantilever_ML (including rotations in null space): "
    os.chdir("cantilever2d.gid")
    sys.path.append(os.getcwd())

    print("Running Trilinos -- Cantilever_ML...")
    success, Msg = benchmarking.MPIParallelRunBenchmark("cantilever_ML.py", "cantilever.txt")

    if(successful==True):
        Text += "OK\n"
        print("Trilinos Cantilever_ML example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("cantilever_ML.py example FAILED")

    os.chdir("..")

    #
    #
    #
    print("resume of all of the examples for the fluid application :")
    print(Text)
    return Text

if __name__ == '__main__':
    Run()
