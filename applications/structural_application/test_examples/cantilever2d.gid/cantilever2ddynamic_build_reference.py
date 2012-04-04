import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

benchmarking.BuildReferenceData("cantilever2ddynamic_benchmarking.py", "cantilever2ddynamic.txt")
