from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
    Msg = ""
    Text = "===== Solid Mechanics Application =====\n"

    print("---start solid mechanics application tests---")

    #
    # SOLID 2D elements tests
   
    # Thin Plate Axial Load
    Text += "2D Small Displacements Plane Stress Axial Load test: "
    os.chdir("1_Thin_plate_Axial_load.gid")
    sys.path.append(os.getcwd())

    print("running the 2D Small Displacements Plane Stress Axial Load test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")

    if(successful==True):
        Text += "OK\n"
        print("2D Small Displacements Plane Stress Axial Load test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("2D Small Displacements Plane Stress Axial Load test FAILED")

    os.chdir("..")
 

    # Thin Plate Self Weight
    Text += "2D Large Displacements Plane Strain Self Weight test: "
    os.chdir("2_Thin_plate_Selfweight.gid")
    sys.path.append(os.getcwd())

    print("running the 2D Large Displacements Plane Strain Self Weight test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")

    if(successful==True):
        Text += "OK\n"
        print("2D Large Displacements Plane Strain Self Weight test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("2D Large Displacements Plane Strain Self Weight test FAILED")

    os.chdir("..")

    # Plane Strain
    Text += "2D Isotropic Damage Plane Stress SimoJu test: "
    os.chdir("IsotropicDamageSimoJuPlaneStress_FourPointShearTest.gid")
    sys.path.append(os.getcwd())

    print("running the 2D Isotropic Damage Plane Stress SimoJu test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")

    if(successful==True):
        Text += "OK\n"
        print("2D Isotropic Damage Plane Stress SimoJu test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("2D Isotropic Damage Plane Stress SimoJu test FAILED")

    os.chdir("..")


    #
    # SOLID 2D elements tests

    # Bending 3D Beam
    Text += "3D Small Displacements Linear Elastic Point Load test: "
    os.chdir("4_Flexion_3D_Beam.gid")
    sys.path.append(os.getcwd())

    print("running the 3D Small Displacements Linear Elastic Point Load test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")

    if(successful==True):
        Text += "OK\n"
        print("3D Small Displacements Linear Elastic Point Load test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("3D Small Displacements Linear Elastic Point Load test FAILED")

    os.chdir("..")
    
    # Z cantilever
    Text += "3D Large Displacements Hyper Elastic Point Load test: "
    os.chdir("4_Flexion_3D_Beam.gid")
    sys.path.append(os.getcwd())

    print("running the 3D Large Displacements Hyper Elastic Point Load test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")

    if(successful==True):
        Text += "OK\n"
        print("3D Large Displacements Hyper Elastic Point Load test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("3D Large Displacements Hyper Elastic Point Load test FAILED")

    os.chdir("..")

    #
    #
    print("---end solid mechanics application tests---")
    print("resume of all of the examples for the SolidMechanicsApplication :")
    print(Text)
    return Text
    

if __name__ == '__main__':
    Run()

