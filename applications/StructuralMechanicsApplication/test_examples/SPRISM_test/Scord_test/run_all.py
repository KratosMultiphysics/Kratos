from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking

def Run():
    Msg = ""
    Text = "===== Structural Mechanics Application - Sprism3D6N Scord test =====\n"

    # Sprism3D6N Scord 2 test:
    Text += "Sprism3D6N Scord 2 test: "
    os.chdir("Sprism3D6N_Scord_2_test.gid")
    sys.path.append(os.getcwd())
    print("running the Sprism3D6N Scord 2 test benchmark...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Sprism3D6N Scord 2 test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Sprism3D6N Scord 2 test  FAILED")
    os.chdir("..")
    
    # Sprism3D6N Scord 4 test:
    Text += "Sprism3D6N Scord 4 test: "
    os.chdir("Sprism3D6N_Scord_4_test.gid")
    sys.path.append(os.getcwd())
    print("running the Sprism3D6N Scord 4 test benchmark...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Sprism3D6N Scord 4 test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Sprism3D6N Scord 4 test  FAILED")
    os.chdir("..")
    
    # Sprism3D6N Scord 8 test:
    Text += "Sprism3D6N Scord 8 test: "
    os.chdir("Sprism3D6N_Scord_8_test.gid")
    sys.path.append(os.getcwd())
    print("running the Sprism3D6N Scord 8 test benchmark...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Sprism3D6N Scord 8 test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Sprism3D6N Scord 8 test  FAILED")
    os.chdir("..")
    
    # Sprism3D6N Scord 16 test:
    Text += "Sprism3D6N Scord 16 test: "
    os.chdir("Sprism3D6N_Scord_16_test.gid")
    sys.path.append(os.getcwd())
    print("running the Sprism3D6N Scord 16 test benchmark...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Sprism3D6N Scord 16 test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Sprism3D6N Scord 16 test  FAILED")
    os.chdir("..")
    
    # Sprism3D6N Scord 32 test:
    Text += "Sprism3D6N Scord 32 test: "
    os.chdir("Sprism3D6N_Scord_32_test.gid")
    sys.path.append(os.getcwd())
    print("running the Sprism3D6N Scord 32 test benchmark...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Sprism3D6N Scord 32 test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Sprism3D6N Scord 32 test  FAILED")
    os.chdir("..")
    
    #
    print("---end forming application tests---")
    print("resume of all of the examples for the Structural Mechanics Application - Sprism3D6N Scord test  :")
    print(Text)
    return Text
    

if __name__ == '__main__':
    Run()

 
