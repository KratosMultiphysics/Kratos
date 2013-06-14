import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for rotating_ball_no_tangent_damp.py..."
benchmarking.BuildReferenceData("rotating_ball_no_tangent_damp.py", "rotating_ball_no_tangent_damp_ref.txt")
