import benchmarking

def Run():
	print "Running two_balls_no_damp test..."
	return benchmarking.RunBenchmark("script.py", "two_balls_no_damp_ref.txt")
