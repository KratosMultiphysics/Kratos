from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
    Msg = ""
    Text = "========== DEM ==========\n"
    Text = "      No benchmarks\n"
    ##

    #Text += "DEM element no damp test: "
    #os.chdir("two_balls_no_damp.gid")
    #sys.path.append(os.getcwd())

    #import update_problem
    #print("Problem updated...")

    #import update_script
    #print("Python script updated...")

    #print("running the DEM two_balls_no_damp test...")
    #successful,Msg = benchmarking.RunBenchmark("two_balls_no_damp.py", "two_balls_no_damp_ref.txt")

    #if(successful==True):
        #Text += "OK\n"
        #print("two_balls_no_damp test successful")
    #else:
        #Text += "FAILED\n"
        #Text += Msg
        #Text += "\n\n"
        #print("two_balls_no_damp example test FAILED")

    #os.chdir("..")
   

    # Add other examples here

    #
    print("resume of all of the examples for the DEM application :")
    print(Text)
    return Text

if __name__ == '__main__':
    Run()
