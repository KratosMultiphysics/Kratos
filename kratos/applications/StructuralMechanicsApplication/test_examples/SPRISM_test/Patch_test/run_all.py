from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking

def Run():
    Msg = ""
    Text = "===== Structural Mechanics Application - Sprism3D6N Patch test =====\n"
 
    # Sprism3D6N Membrane patch test:
    Text += "Sprism3D6N Membrane patch test: "
    os.chdir("Sprism3D6N_membrane_patch_test.gid")
    sys.path.append(os.getcwd())
    print("running the Sprism3D6N Membrane patch test benchmark...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Sprism3D6N Membrane patch test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Sprism3D6N Membrane patch test  FAILED")
    os.chdir("..")
    
    # Sprism3D6N Bendig patch test:
    Text += "Sprism3D6N Bending patch test: "
    os.chdir("Sprism3D6N_bending_patch_test.gid")
    sys.path.append(os.getcwd())
    print("running the Sprism3D6N Bending patch test benchmark...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Sprism3D6N Bending patch test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Sprism3D6N Bending patch test FAILED")
    os.chdir("..")
    
    #
    print("---end forming application tests---")
    print("resume of all of the examples for the Structural Mechanics Application - Sprism3D6N Patch test  :")
    print(Text)
    return Text
    

if __name__ == '__main__':
    Run()

 
