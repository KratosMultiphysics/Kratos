import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for testConvection.py..."
benchmarking.BuildReferenceData("test_pureconvectionsolver_benchmarking.py", "test_pureconvectionsolver_benchmarking_ref.txt")
