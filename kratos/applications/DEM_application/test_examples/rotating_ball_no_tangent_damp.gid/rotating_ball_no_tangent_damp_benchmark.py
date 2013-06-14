import benchmarking

def Run():
	print "Running rotating_ball_no_tangent_damp test..."
	return benchmarking.RunBenchmark("rotating_ball_no_tangent_damp.py", "rotating_ball_no_tangent_damp_ref.txt")
