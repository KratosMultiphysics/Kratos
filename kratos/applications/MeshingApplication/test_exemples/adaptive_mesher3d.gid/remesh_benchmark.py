kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking

import sys
sys.path.append(kratos_benchmarking_path)

import benchmarking

def Run():        
	print "Running remesh.py..."
	return benchmarking.RunBenchmark("remesh.py", "adaptive_mesher3d_ref.txt")

