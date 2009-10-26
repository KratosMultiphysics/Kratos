import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for Moving_Circle..."
benchmarking.BuildReferenceData("ProjectionTest_2D_benchmarking.py", "ProjectionTest_2D_benchmarking_ref.txt")
