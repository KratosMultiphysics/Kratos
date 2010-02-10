import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "running the benchmark for mass_conservation test..."
Msg = benchmarking.RunBenchmark("run_benchmark.py", "benchmark_reference_solution.txt")

if (Msg == True):
    print "mass_conservation test example succesful"
else:
    print "mass_conservation test example FAILED"
