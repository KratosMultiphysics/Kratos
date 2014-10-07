from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
kratos_benchmarking_path = '../../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking

Text = " "

#
print("verifying  Couette2dNonNewt.py...")
successful,Msg = benchmarking.RunBenchmark("Couette2dNonNewt.py", "Couette2dNonNewt_ref.txt")

if(successful==True):
    Text += "OK\n"
    print("Couette2dNonNewt example SUCCESFUL")
else:
    Text += "FAILED\n"
    Text += Msg
    Text += "\n\n"
    print("Couette2dNonNewt example FAILED")

print(Text)
