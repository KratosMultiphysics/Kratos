import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for Patch_Test_Total_Lagrangian_4N_benchmarking.py..."
benchmarking.BuildReferenceData("Patch_Test_Total_Lagrangian_4N_benchmarking.py", "Patch_Test_Total_Lagrangian_4N.txt")
