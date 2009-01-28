import benchmarking

def Run():
	print "Running dam2d.py..."
	return benchmarking.RunBenchmark("dam2d.py", "dam2d_ref.txt")

