import benchmarking

def Run():
	print "Running rotating_ball_rolling_friction test..."
	return benchmarking.RunBenchmark("rotating_ball_rolling_friction.py", "rotating_ball_rolling_friction_ref.txt")
