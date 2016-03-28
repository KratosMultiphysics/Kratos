from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking

def Run():
    Msg = ""
    Text = "===== Structural Mechanics Application - Sprism3D6N Sphere test =====\n"
 
    # Sprism3D6N Sphere 4 test:
    Text += "Sprism3D6N Sphere 4 test: "
    os.chdir("Sprism3D6N_Sphere_4_test.gid")
    sys.path.append(os.getcwd())
    print("running the Sprism3D6N Sphere 4 test benchmark...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Sprism3D6N Sphere 4 test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Sprism3D6N Sphere 4 test  FAILED")
    os.chdir("..")
    
    # Sprism3D6N Sphere 8 test:
    Text += "Sprism3D6N Sphere 8 test: "
    os.chdir("Sprism3D6N_Sphere_8_test.gid")
    sys.path.append(os.getcwd())
    print("running the Sprism3D6N Sphere 8 test benchmark...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Sprism3D6N Sphere 8 test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Sprism3D6N Sphere 8 test  FAILED")
    os.chdir("..")
    
    # Sprism3D6N Sphere 16 test:
    Text += "Sprism3D6N Sphere 16 test: "
    os.chdir("Sprism3D6N_Sphere_16_test.gid")
    sys.path.append(os.getcwd())
    print("running the Sprism3D6N Sphere 16 test benchmark...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Sprism3D6N Sphere 16 test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Sprism3D6N Sphere 16 test  FAILED")
    os.chdir("..")
    
    # Sprism3D6N Sphere 24 test:
    Text += "Sprism3D6N Sphere 24 test: "
    os.chdir("Sprism3D6N_Sphere_24_test.gid")
    sys.path.append(os.getcwd())
    print("running the Sprism3D6N Sphere 24 test benchmark...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Sprism3D6N Sphere 24 test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Sprism3D6N Sphere 24 test  FAILED")
    os.chdir("..")
    
    # Sprism3D6N Sphere 32 test:
    Text += "Sprism3D6N Sphere 32 test: "
    os.chdir("Sprism3D6N_Sphere_32_test.gid")
    sys.path.append(os.getcwd())
    print("running the Sprism3D6N Sphere 32 test benchmark...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Sprism3D6N Sphere 32 test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Sprism3D6N Sphere 32 test  FAILED")
    os.chdir("..")
    
    #
    print("---end forming application tests---")
    print("resume of all of the examples for the Structural Mechanics Application - Sprism3D6N Sphere test  :")
    print(Text)
    return Text
    

if __name__ == '__main__':
    Run()

 
