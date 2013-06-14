import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for two_balls_normal_damp.py..."
benchmarking.BuildReferenceData("two_balls_normal_damp.py", "two_balls_normal_damp_ref.txt")
