import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for QuietWater_LevelSet.py..."
benchmarking.BuildReferenceData("split_level_set_QuietWater.py", "QuietWater_LS_ref.txt")
