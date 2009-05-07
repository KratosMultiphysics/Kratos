import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for cil_gls.py..."
benchmarking.BuildReferenceData("cil_gls.py", "cil_gls_ref.txt")
