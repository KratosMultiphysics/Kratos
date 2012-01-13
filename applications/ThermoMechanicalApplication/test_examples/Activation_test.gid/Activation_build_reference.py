import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for Activation_test..."
benchmarking.BuildReferenceData("run_benchmark.py", "Activation_reference_solution.txt")
