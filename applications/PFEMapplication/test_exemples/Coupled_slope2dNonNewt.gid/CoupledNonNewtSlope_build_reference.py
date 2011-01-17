import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for CoupledNonNewtSlope.py..."
benchmarking.BuildReferenceData("CoupledNonNewtSlope.py", "CoupledNonNewtSlope_ref.txt")
