from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
    Msg = ""
    Text = "== Solid Mechanics ==========\n"

    #
    # VMS2D element test

    Text += "scordelis low roof test: "
    os.chdir("scordelis.gid")
    sys.path.append(os.getcwd())

    print("running the scordelis low roof benchmark test...")
    Msg = benchmarking.RunBenchmark("run_test.py", "min_displacements.txt")

    if (Msg):
        Text += "OK\n"
        print("scordelis low roof test succesful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("scordelis low roof test FAILED")

    os.chdir("..")

    #
    
    #
    # SHELL elements tests
    
    Text += "Shell Q4 Thick Bending RollUp test: "
    os.chdir("Shell_Q4_Thick__BendingRollUp.gid")
    sys.path.append(os.getcwd())
    print("running the Shell Q4 Thick Bending RollUp benchmark test...")
    Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if (Msg):
        Text += "OK\n"
        print("Shell Q4 Thick Bending RollUp test succesful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Shell Q4 Thick Bending RollUp test FAILED")
    os.chdir("..")
    
    Text += "Shell Q4 Thick Drilling RollUp test: "
    os.chdir("Shell_Q4_Thick__DrillingRollUp.gid.gid")
    sys.path.append(os.getcwd())
    print("running the Shell Q4 Thick Drilling RollUp benchmark test...")
    Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if (Msg):
        Text += "OK\n"
        print("Shell Q4 Thick Drilling RollUp test succesful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Shell Q4 Thick Drilling RollUp test FAILED")
    os.chdir("..")
    
    Text += "Shell T3 Thin Bending RollUp test: "
    os.chdir("Shell_T3_Thin__BendingRollUp.gid")
    sys.path.append(os.getcwd())
    print("running the Shell T3 Thin Bending RollUp benchmark test...")
    Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if (Msg):
        Text += "OK\n"
        print("Shell T3 Thin Bending RollUp test succesful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Shell T3 Thin Bending RollUp test FAILED")
    os.chdir("..")
    
    Text += "Shell T3 Thin Drilling RollUp test: "
    os.chdir("Shell_T3_Thin__DrillingRollUp.gid.gid")
    sys.path.append(os.getcwd())
    print("running the Shell T3 Thin Drilling RollUp benchmark test...")
    Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if (Msg):
        Text += "OK\n"
        print("Shell T3 Thin Drilling RollUp test succesful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Shell T3 Thin Drilling RollUp test FAILED")
    os.chdir("..")
    
    #
    #
    
    print("resume of all of the examples for the Solid Mechanics application :")
    print(Text)
    return Text
    

if __name__ == '__main__':
    Run()
