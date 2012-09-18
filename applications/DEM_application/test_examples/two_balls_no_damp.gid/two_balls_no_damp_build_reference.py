import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for script.py..."
benchmarking.BuildReferenceData("two_balls_no_damp.py", "two_balls_no_damp_ref.txt")
