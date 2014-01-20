from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
kratos_benchmarking_path = '../../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking

Text = " "

#
print("verifying  CoupledNonNewtSlope.py...")
Msg = benchmarking.RunBenchmark("CoupledNonNewtSlope.py", "CoupledNonNewtSlope_ref.txt")

if (Msg):
    Text += "OK\n"
    print("CoupledNonNewtSlope example SUCCESFUL")
else:
    Text += "FAILED\n"
    Text += Msg
    Text += "\n\n"
    print("CoupledNonNewtSlope example FAILED")

print(Text)
