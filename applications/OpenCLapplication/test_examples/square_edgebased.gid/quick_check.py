import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "running reference data for edgebased_PureConvection.py..."
benchmarking.RunBenchmark("edgebased_PureConvection.py", "test_pureconvectionsolver_benchmarking_ref.txt")
