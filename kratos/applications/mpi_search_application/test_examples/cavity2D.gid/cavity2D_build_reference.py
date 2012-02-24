import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for test_fractstep_cont_laplacian.py..."
benchmarking.BuildReferenceData("test_fractstep_cont_laplacian.py", "fractstep_cont_laplacian_benchmarking_ref.txt")

print "Building reference data for test_fractstep_discrete_laplacian.py..."
benchmarking.BuildReferenceData("test_fractstep_discrete_laplacian.py", "fractstep_discrete_laplacian_benchmarking_ref.txt")
