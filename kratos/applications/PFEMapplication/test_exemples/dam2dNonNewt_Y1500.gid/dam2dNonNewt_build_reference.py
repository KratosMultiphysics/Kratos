import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for dam2dNonNewt.py..."
benchmarking.BuildReferenceData("dam2dNonNewt.py", "dam2dNonNewt_ref.txt")
