from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
    Msg = ""
    Text = "========== Meshing Aplication ==========\n"

    #

    Text += "adaptive_mesher2d: "
    os.chdir("adaptive_mesher2d.gid")
    sys.path.append(os.getcwd())

    print("Running Adaptive Mesher 2d benchmark...")
    successful,Msg = benchmarking.RunBenchmark("remesh.py", "adaptive_mesher2d_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("Adaptive_mesher2d benchmarking example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Adaptive_mesher2d benchmarking example FAILED")

    os.chdir("..")

    #

    Text += "adaptive_mesher3d: "
    os.chdir("adaptive_mesher3d.gid")
    sys.path.append(os.getcwd())

    print("Running Adaptive Mesher 3d benchmark...")
    successful,Msg = benchmarking.RunBenchmark("remesh.py", "adaptive_mesher3d_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("Adaptive_mesher2d benchmarking example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Adaptive_mesher2d benchmarking example FAILED")

    os.chdir("..")

    #
    Text += "Mapping_2d: "
    os.chdir("Mapping_2d.gid")
    sys.path.append(os.getcwd())

    print("Running Mapping_2d benchmark...")
    successful,Msg = benchmarking.RunBenchmark("ProjectionTest_2D_script.py", "ProjectionTest_2D_benchmarking_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("Mapping_2d benchmarking example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Mapping_2d benchmarking example FAILED")

    os.chdir("..")

    #
    Text += "Mapping_3d: "
    os.chdir("Mapping_3d.gid")
    sys.path.append(os.getcwd())

    print("Running Mapping_3d benchmark...")
    successful,Msg = benchmarking.RunBenchmark("ProjectionTest_3D_script.py", "ProjectionTest_3D_benchmarking_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("Mapping_3d benchmarking example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Mapping_3d benchmarking example FAILED")

    os.chdir("..")

    #

    Text += "Mapping_2d_BinBased: "
    os.chdir("Mapping_2d_BinBased.gid")
    sys.path.append(os.getcwd())

    print("Running Mapping_2d_BinBased benchmark...")
    successful,Msg = benchmarking.RunBenchmark("ProjectionTestBinBased_2D_script.py", "ProjectionTestBinBased_2D_benchmarking_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("Mapping_2d_BinBased benchmarking example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Mapping_2d_BinBased benchmarking example FAILED")

    os.chdir("..")

   #

    Text += "Mapping_3d_BinBased: "
    os.chdir("Mapping_3d_BinBased.gid")
    sys.path.append(os.getcwd())

    print("Running Mapping_3d_BinBased benchmark...")
    successful,Msg = benchmarking.RunBenchmark("ProjectionTestBinBased_3D_script.py", "ProjectionTestBinBased_3D_benchmarking_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("Mapping_3d_BinBased benchmarking example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Mapping_3d_BinBased benchmarking example FAILED")

    os.chdir("..")

   #
    Text += "connectivity_preserve_modeler benchmark: "
    os.chdir("test_connectivity_preserve_modeler")
    sys.path.append(os.getcwd())
    print("verifying  connectivity_preserve_modeler benchmark...")
    successful,Msg = benchmarking.RunBenchmark("do_test.py", "connectivity_preserve_benchmarking_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("connectivity_preserve_modeler benchmark successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("connectivity_preserve_modeler benchmark FAILED")

    os.chdir("..")

    #
    #

    print("resume of all of the examples for the meshing application :")
    print(Text)
    return Text

if __name__ == '__main__':
    Run()
