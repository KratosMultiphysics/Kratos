import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for naca3D test..."
benchmarking.BuildReferenceData("run_benchmark.py", "benchmark_reference_solution.txt")
