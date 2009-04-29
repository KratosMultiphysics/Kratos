import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for balken.gid..."
benchmarking.BuildReferenceData("balken_benchmarking.py", "balken_ref.txt")
