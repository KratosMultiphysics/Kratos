import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for Mapping 2d..."
benchmarking.BuildReferenceData("ProjectionTest_2D_benchmarking.py", "ProjectionTest_2D_benchmarking_ref.txt")
