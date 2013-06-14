import benchmarking

def Run():
	print "Running two_balls_normal_damp test..."
	return benchmarking.RunBenchmark("two_balls_normal_damp.py", "two_balls_normal_damp_ref.txt")
