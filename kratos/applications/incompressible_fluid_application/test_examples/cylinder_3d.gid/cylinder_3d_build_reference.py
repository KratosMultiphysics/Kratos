import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for cylinder_3d.py..."
benchmarking.BuildReferenceData("run_example.py", "cylinder_3d_ref.txt")
