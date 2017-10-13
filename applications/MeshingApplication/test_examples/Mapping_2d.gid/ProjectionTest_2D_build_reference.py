from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
kratos_benchmarking_path = '../../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking

print("Building reference data for Mapping 2d...")
benchmarking.BuildReferenceData("ProjectionTest_2D_script.py", "ProjectionTest_2D_benchmarking_ref.txt")
