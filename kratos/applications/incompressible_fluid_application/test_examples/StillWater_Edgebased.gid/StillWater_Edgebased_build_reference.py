import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for StillWater_Edgebased.py..."
benchmarking.BuildReferenceData("StillWater_Edgebased_script.py", "StillWater_Edgebased_ref.txt")
