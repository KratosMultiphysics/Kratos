# Import Kratos core and apps
import KratosMultiphysics as KM

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
import os, csv

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

model = KM.Model()

# Create optimizer and perform optimization
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model)
optimizer.Optimize()

# =======================================================================================================
# Test results and clean directory
# =======================================================================================================
original_directory = os.getcwd()
output_directory = parameters["optimization_settings"]["output"]["output_directory"].GetString()
optimization_model_part_name = parameters["optimization_settings"]["model_settings"]["model_part_name"].GetString()
optimization_log_filename = parameters["optimization_settings"]["output"]["optimization_log_filename"].GetString() + ".csv"

# Testing by
# 1) using the "json_output_process" & "json_check_process" within the structural analysis
# 2) additionally checking some process output
os.chdir(output_directory)

with open(optimization_log_filename, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    last_line = None
    for line in reader:
        if not line:
            continue
        else:
            last_line = line

    resulting_iteration = float(last_line[0].strip())
    resulting_abs_improvement = float(last_line[2].strip())

    # Check against specifications
    TestCase().assertEqual(resulting_iteration, 7)
    TestCase().assertAlmostEqual(resulting_abs_improvement, -98.843, 3)

os.chdir(original_directory)

# Cleaning
kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
kratos_utilities.DeleteDirectoryIfExisting(output_directory)
kratos_utilities.DeleteFileIfExisting(os.path.basename(original_directory)+".post.lst")
kratos_utilities.DeleteFileIfExisting(optimization_model_part_name+".time")

# =======================================================================================================