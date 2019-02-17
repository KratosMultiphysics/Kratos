from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
    Msg = ""
    Text = "========== Structural Aplication ==========\n"

    #

    Text += "Patch_Test_Total_Lagrangian_3D_8N: "
    os.chdir("Patch_Test_Total_Lagrangian_3D_8N.gid")
    sys.path.append(os.getcwd())

    print("Running Patch_Test_Total_Lagrangian_3D_8N...")
    successful,Msg = benchmarking.RunBenchmark("Patch_Test_Total_Lagrangian_3D_8N_benchmarking.py", "Patch_Test_Total_Lagrangian_3D_8N.txt")

    if(successful==True):
        Text += "OK\n"
        print("Patch_Test_Total_Lagrangian_3D_8N_benchmarking example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Patch_Test_Total_Lagrangian_3D_8N_benchmarking.py example FAILED")

    os.chdir("..")

    #
    Text += "Patch_Test_Total_Lagrangian_4N: "
    os.chdir("Patch_Test_Total_Lagrangian_4N.gid")
    sys.path.append(os.getcwd())

    print("Running Patch_Test_Total_Lagrangian_4N...")
    successful,Msg = benchmarking.RunBenchmark("Patch_Test_Total_Lagrangian_4N_benchmarking.py", "Patch_Test_Total_Lagrangian_4N.txt")

    if(successful==True):
        Text += "OK\n"
        print("Patch_Test_Total_Lagrangian_4N_benchmarking example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Patch_Test_Total_Lagrangian_4N_benchmarking.py example FAILED")

    os.chdir("..")

   #
    Text += "cantilever2dDynamic:"
    os.chdir("cantilever2d.gid")
    sys.path.append(os.getcwd())

    print("Running cantilever2ddynamic")
    successful,Msg = benchmarking.RunBenchmark("cantilever2ddynamic_benchmarking.py", "cantilever2ddynamic.txt")

    if(successful==True):
        Text += "OK\n"
        print("cantilever2dynamic successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("cantilever2dynamic example FAILED")

    os.chdir("..")

    #

    Text += "cantilever2dstatic using MKLPardisoSolver: "
    os.chdir("cantilever2d.gid")
    sys.path.append(os.getcwd())

    print("Running cantilever2dstatic using MKLPardisoSolver")
    successful,Msg = benchmarking.RunBenchmark("cantilever2dstatic_Mkl_benchmarking.py", "cantilever2dstatic_superlu.txt")

    if(successful==True):
        Text += "OK\n"
        print("cantilever2dstatic using MKLPardisoSolver successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Expected fail because the  program should be compiled with INTEL compliler.")

    os.chdir("..")

    #

    Text += "cantilever2dstatic_superlu: "
    os.chdir("cantilever2d.gid")
    sys.path.append(os.getcwd())

    print("Running cantilever2dstatic witdh SuperLUSolver")
    successful,Msg = benchmarking.RunBenchmark("cantilever2dstatic_superlu_benchmarking.py", "cantilever2dstatic_superlu.txt")

    if(successful==True):
        Text += "OK\n"
        print("cantilever2dstatic witdh SuperLUSolver example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("cantilever2dstatic witdh SuperLUSolver example FAILED")

    os.chdir("..")

    #
    Text += "cantilever2dstatic: "
    os.chdir("cantilever2d.gid")
    sys.path.append(os.getcwd())

    print("Running cantilever2dstatic with classical solver")
    successful,Msg = benchmarking.RunBenchmark("cantilever2dstatic_benchmarking.py", "cantilever2dstatic_superlu.txt")

    if(successful==True):
        Text += "OK\n"
        print("cantilever2dstatic with classical solver example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("cantilever2dstatic with classical solver example FAILED")

    os.chdir("..")
    #

    Text += "cantilever2dstatic_parallel:"
    os.chdir("cantilever2d.gid")
    sys.path.append(os.getcwd())

    print("Running cantilever2dstatic_parallel using ParallelSkylineLUFactorizationSolver")
    successful,Msg = benchmarking.RunBenchmark("cantilever2dstatic_parallel_benchmarking.py", "cantilever2dstatic_superlu.txt")

    if(successful==True):
        Text += "OK\n"
        print("cantilever2dstatic_parallel using ParallelSkylineLUFactorizationSolver  example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Expected fail because the  program should be compiled with INTEL compliler.")

    os.chdir("..")

    #
    Text += "cantilever2dstatic using ParallelMKLPardisoSolver: "
    os.chdir("cantilever2d.gid")
    sys.path.append(os.getcwd())

    print("Running cantilever2dstatic using ParallelMKLPardisoSolver")
    successful,Msg = benchmarking.RunBenchmark("cantilever2dstatic_parallel_pardiso_benchmarking.py", "cantilever2dstatic_superlu.txt")

    if(successful==True):
        Text += "OK\n"
        print("cantilever2dstatic using ParallelMKLPardisoSolver  example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Expected fail because the  program should be compiled with INTEL compliler.")

    os.chdir("..")

    #

    Text += "cantilever3dstatic: "
    os.chdir("cantilever3d.gid")
    sys.path.append(os.getcwd())

    print("Running cantilever3dstatic")
    successful,Msg = benchmarking.RunBenchmark("cantilever3dstatic_superlu_benchmarking.py", "cantilever3dstatic.txt")

    if(successful==True):
        Text += "OK\n"
        print("cantilever3dstatic example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("cantilever3dstatic example FAILED")

    os.chdir("..")

    #
    Text += "cantilever3ddynamic: "
    os.chdir("cantilever3d.gid")
    sys.path.append(os.getcwd())

    print("Running cantilever3dynamic")
    successful,Msg = benchmarking.RunBenchmark("cantilever3ddynamic_benchmarking.py", "cantilever3ddynamic.txt")

    if(successful==True):
        Text += "OK\n"
        print("cantilever3ddynamic example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("cantilever2dynamic  example FAILED")

    os.chdir("..")

    #

    Text += "balken contact benchmark: "
    os.chdir("balken.gid")
    sys.path.append(os.getcwd())

    print("Running balken contact benchmark...")
    successful,Msg = benchmarking.RunBenchmark("balken_benchmarking.py", "balken_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("balken contact benchmark example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("balken contact benchmark example FAILED")

    os.chdir("..")

    #
    Text += "arc length desplacement benchmark: "
    os.chdir("arc_length_des.gid")
    sys.path.append(os.getcwd())

    print("Running arc length desplacement  benchmark...")
    successful,Msg = benchmarking.RunBenchmark("arc_length_des_benchmarking.py", "arc_length_des_benchmarking.txt")

    if(successful==True):
        Text += "OK\n"
        print("arc_length_des benchmark example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("arc_length_des benchmark example FAILED")

    os.chdir("..")

    #
    Text += "arc length benchmark: "
    os.chdir("arc_length.gid")
    sys.path.append(os.getcwd())

    print("Running arc length  benchmark...")
    successful,Msg = benchmarking.RunBenchmark("arc_length_benchmarking.py", "arc_length_benchmarking.txt")

    if(successful==True):
        Text += "OK\n"
        print("arc_length benchmark example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("arc_length benchmark example FAILED")

    os.chdir("..")

    #
    Text += "Pendulo Kratos Length benchmark: "
    os.chdir("Pendulo_Kratos.gid")
    sys.path.append(os.getcwd())

    print("Pendulo Kratos Length benchmark...")
    successful,Msg = benchmarking.RunBenchmark("Pendulo_Kratos_benchmarking.py", "Pendulo_Kratos_benchmarking.txt")

    if(successful==True):
        Text += "OK\n"
        print("Pendulo Kratos benchmark example successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Pendulo Kratos benchmark example FAILED")

    os.chdir("..")

    #

    Text += "PlasticitJ2 with EBST collection benchmarks: "
    os.chdir("plasticJ2.gid")
    sys.path.append(os.getcwd())

    print("PlasticitJ2 with EBST TENSION benchmarks...")
    successful,Msg = benchmarking.RunBenchmark("tension.py", "tension_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("PlasticitJ2 with EBST TENSION benchmarks successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("PlasticitJ2 with EBST TENSION benchmarks FAILED")

    print("PlasticitJ2 with EBST TORSION benchmarks...")
    successful,Msg = benchmarking.RunBenchmark("torsion.py", "torsion_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("PlasticitJ2 with EBST TORSION benchmarks successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("PlasticitJ2 with EBST TORSION benchmarks FAILED")

    print("PlasticitJ2 with EBST VERTICAL benchmarks...")
    successful,Msg = benchmarking.RunBenchmark("vertical.py", "vertical_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("PlasticitJ2 with EBST VERTICAL benchmarks successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("PlasticitJ2 with EBST VERTICAL benchmarks FAILED")

    print("PlasticitJ2 with EBST FORCE benchmarks...")
    successful,Msg = benchmarking.RunBenchmark("force.py", "force_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("PlasticitJ2 with EBST FORCE benchmarks successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("PlasticitJ2 with EBST FORCE benchmarks FAILED")

    os.chdir("..")

    #

    print("Resume of all of the examples for the structural application :")
    print(Text)
    return Text

if __name__ == '__main__':
    Run()
