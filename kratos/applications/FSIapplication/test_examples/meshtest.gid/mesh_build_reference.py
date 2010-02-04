import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for mesh.py..."
benchmarking.BuildReferenceData("mesh.py", "mesh_ref.txt")
