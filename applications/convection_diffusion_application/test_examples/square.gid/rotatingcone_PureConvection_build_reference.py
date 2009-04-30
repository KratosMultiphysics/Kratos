import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for rotatingcone_PureConvection.py..."
benchmarking.BuildReferenceData("rotatingcone_PureConvectionBenchmarking.py", "rotatingcone_PureConvection_ref.txt")
