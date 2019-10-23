# Import Kratos core and apps
import KratosMultiphysics as KM

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
import csv, os

# Read parameters
with open("parameters.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

model = KM.Model()

# Create optimizer and perform optimization
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model)
optimizer.Optimize()

# =======================================================================================================
# Test results and clean directory
# =======================================================================================================
output_directory = parameters["optimization_settings"]["output"]["output_directory"].GetString()
optimization_log_filename = parameters["optimization_settings"]["output"]["optimization_log_filename"].GetString() + ".csv"
optimization_model_part_name = parameters["optimization_settings"]["model_settings"]["model_part_name"].GetString()

# Testing
original_directory = os.getcwd()
os.chdir(output_directory)

with open(optimization_log_filename, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    last_line = None
    for line in reader:
        if not line:
            continue
        else:
            last_line = line

    resulting_optimization_iterations = int(last_line[0].strip())
    resulting_improvement = float(last_line[2].strip())

    # Check against specifications
    TestCase().assertEqual(resulting_optimization_iterations, 16)
    TestCase().assertAlmostEqual(resulting_improvement, -17.553770, 4)

os.chdir(original_directory)

# Cleaning
kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
kratos_utilities.DeleteDirectoryIfExisting(output_directory)
kratos_utilities.DeleteFileIfExisting(os.path.basename(original_directory)+".post.lst")
kratos_utilities.DeleteFileIfExisting(optimization_model_part_name+".time")

# =======================================================================================================