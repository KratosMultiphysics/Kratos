import benchmarking

def Run():
	print "Running slip condition test..."
	return benchmarking.RunBenchmark("slip_test.py", "slip_test_ref.txt")
