import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "running the benchmark for naca3D test..."
Msg = benchmarking.RunBenchmark("run_benchmark.py", "benchmark_reference_solution.txt")

if (Msg == True):
    print "naca3D test example succesful"
else:
    print "naca3D test example FAILED"
