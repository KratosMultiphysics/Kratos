import sys

kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)

import benchmarking

print "Building reference data for slip condition test..."
benchmarking.BuildReferenceData("slip_test.py", "slip_test_ref.txt")
