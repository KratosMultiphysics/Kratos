import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for Mapping_3d..."
benchmarking.BuildReferenceData("ProjectionTestBinBased_3D_script.py", "ProjectionTestBinBased_3D_benchmarking_ref.txt")
