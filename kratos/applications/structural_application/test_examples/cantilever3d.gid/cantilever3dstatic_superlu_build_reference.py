import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

benchmarking.BuildReferenceData("cantilever3dstatic_superlu_benchmarking.py", "cantilever3dstatic.txt")
