import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for Couette2dNonNewt.py..."
benchmarking.BuildReferenceData("Couette2dNonNewt.py", "Couette2dNonNewt_ref.txt")
