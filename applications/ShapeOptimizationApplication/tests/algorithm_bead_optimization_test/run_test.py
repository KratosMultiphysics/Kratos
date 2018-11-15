# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
import csv, os

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters = Parameters(parameter_file.read())

model = Model()

# Create optimizer and perform optimization
import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model)
optimizer.Optimize()

# =======================================================================================================
# Test results and clean directory
# =======================================================================================================
output_directory = parameters["optimization_settings"]["output"]["output_directory"].GetString()
response_log_filename = parameters["optimization_settings"]["output"]["response_log_filename"].GetString() + ".csv"
optimization_model_part_name = parameters["optimization_settings"]["model_settings"]["model_part_name"].GetString()

# Testing
original_directory = os.getcwd()
os.chdir(output_directory)

with open(response_log_filename, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    last_line = None
    for line in reader:
        if not line:
            continue
        else:
            last_line = line

    resulting_lagrange_value = float(last_line[3].strip())
    resulting_objective_value = float(last_line[5].strip())
    resulting_lambda = float(last_line[8].strip())
    resulting_penalty_value = float(last_line[9].strip())
    resulting_penalty_scaling = float(last_line[10].strip())

    # Check against specifications
    TestCase().assertAlmostEqual(resulting_lagrange_value, 4.45863E-04,5)
    TestCase().assertAlmostEqual(resulting_objective_value, 4.45786E-04,5)
    TestCase().assertAlmostEqual(resulting_lambda, 1.87765E-03,5)
    TestCase().assertAlmostEqual(resulting_penalty_value, 2.47040E-05,5)
    TestCase().assertAlmostEqual(resulting_penalty_scaling, 3.89608E-07,5)

os.chdir(original_directory)

# Cleaning
kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
kratos_utilities.DeleteDirectoryIfExisting(output_directory)
kratos_utilities.DeleteFileIfExisting(os.path.basename(original_directory)+".post.lst")
kratos_utilities.DeleteFileIfExisting(optimization_model_part_name+".time")
kratos_utilities.DeleteFileIfExisting(optimization_model_part_name+".post.bin")