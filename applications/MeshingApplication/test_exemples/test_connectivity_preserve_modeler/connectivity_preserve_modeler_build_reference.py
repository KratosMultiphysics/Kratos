import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for testing the connectivity_preserve_modeler..."
benchmarking.BuildReferenceData("do_test.py", "connectivity_preserve_benchmarking_ref.txt")
