import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for cavity2d.py..."
benchmarking.BuildReferenceData("cavity2d.py", "cavity_ref.txt")
