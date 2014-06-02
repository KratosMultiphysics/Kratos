from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
    Msg = ""
    Text = "===== Solid Mechanics Application =====\n"

    
    #
    # SHELL elements tests
   
    Text += "Shell T3 Isotropic  Scordelis  test: "
    os.chdir("Shell_T3_Isotropic_Scordelis.gid")
    sys.path.append(os.getcwd())

    print("---start solid mechanics application tests---")

    print("running the Scordelis Low Roof benchmark test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "min_displacements.txt")

    if(successful==True):
        Text += "OK\n"
        print("Scordelis Low Roof test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Scordelis Low Roof test FAILED")

    os.chdir("..")
 
    # Bending RollUp Q4 test:
    Text += "Shell Q4 Thick Bending  RollUp test: "
    os.chdir("Shell_Q4_Thick__BendingRollUp.gid")
    sys.path.append(os.getcwd())
    print("running the Shell Q4 Thick Bending RollUp benchmark test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Shell Q4 Thick Bending RollUp test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Shell Q4 Thick Bending RollUp test FAILED")
    os.chdir("..")

    # Drilling RollUp Q4 test:    
    Text += "Shell Q4 Thick Drilling RollUp test: "
    os.chdir("Shell_Q4_Thick__DrillingRollUp.gid")
    sys.path.append(os.getcwd())
    print("running the Shell Q4 Thick Drilling RollUp benchmark test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Shell Q4 Thick Drilling RollUp test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Shell Q4 Thick Drilling RollUp test FAILED")
    os.chdir("..")
    
    # Bending RollUp T3 test:
    Text += "Shell T3 Thin  Bending  RollUp test: "
    os.chdir("Shell_T3_Thin__BendingRollUp.gid")
    sys.path.append(os.getcwd())
    print("running the Shell T3 Thin Bending RollUp benchmark test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Shell T3 Thin Bending RollUp test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Shell T3 Thin Bending RollUp test FAILED")
    os.chdir("..")
 
    # Drilling RollUp T3 test:  
    Text += "Shell T3 Thin  Drilling RollUp test: "
    os.chdir("Shell_T3_Thin__DrillingRollUp.gid")
    sys.path.append(os.getcwd())
    print("running the Shell T3 Thin Drilling RollUp benchmark test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if(successful==True):
        Text += "OK\n"
        print("Shell T3 Thin Drilling RollUp test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Shell T3 Thin Drilling RollUp test FAILED")
    os.chdir("..")
    
    #
    #
    print("---end solid mechanics application tests---")
    print("resume of all of the examples for the SolidMechanicsApplication :")
    print(Text)
    return Text
    

if __name__ == '__main__':
    Run()

