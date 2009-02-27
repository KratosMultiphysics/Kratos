import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

print "Building reference data for dam2d.py..."
benchmarking.BuildReferenceData("column.py", "column_ref.txt")
