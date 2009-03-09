import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for remesh.py..."
benchmarking.BuildReferenceData("remesh.py", "adapt_remesh3d_ref.txt")
