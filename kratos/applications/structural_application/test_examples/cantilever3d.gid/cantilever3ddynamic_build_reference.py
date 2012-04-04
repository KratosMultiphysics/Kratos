import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

benchmarking.BuildReferenceData("cantilever3ddynamic_benchmarking.py", "cantilever3ddynamic.txt")
