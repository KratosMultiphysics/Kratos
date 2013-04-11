import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for fractional step element: cavity flow test..."
benchmarking.BuildReferenceData("fs_benchmark.py", "fs_benchmark_ref.txt")

