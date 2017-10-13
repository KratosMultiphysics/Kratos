from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
kratos_benchmarking_path = '../../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking

print("Building reference data for cavity_monolithic_3d...")
benchmarking.BuildReferenceData("script.py", "cavity_monolithic_3d_ref.txt")
