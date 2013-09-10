import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for scordelis low roof"
benchmarking.BuildReferenceData("run_test.py", "min_displacements.txt")
