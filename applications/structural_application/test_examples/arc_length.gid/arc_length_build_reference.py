import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for arc_length_des_benchmarking.py..."
benchmarking.BuildReferenceData("arc_length_benchmarking.py",  "arc_length_benchmarking.txt")
