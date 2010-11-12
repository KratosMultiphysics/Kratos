import benchmarking

def Run():
	print "Running VMS2D test..."
	return benchmarking.RunBenchmark("script_elemtest.py", "vms2d_test_ref.txt")
