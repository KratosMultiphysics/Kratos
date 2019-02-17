from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
kratos_benchmarking_path = '../../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking

print("Building reference data for tension.py...")
benchmarking.BuildReferenceData("tension.py", "tension_ref.txt")

print("Building reference data for torsion.py...")
benchmarking.BuildReferenceData("torsion.py", "torsion_ref.txt")

print("Building reference data for vertical.py...")
benchmarking.BuildReferenceData("vertical.py", "vertical_ref.txt")

print("Building reference data for force.py...")
benchmarking.BuildReferenceData("force.py", "force_ref.txt")
# print "Building reference data for test_fractstep_discrete_laplacian.py..."
# benchmarking.BuildReferenceData("test_fractstep_discrete_laplacian.py", "fractstep_discrete_laplacian_benchmarking_ref.txt")
