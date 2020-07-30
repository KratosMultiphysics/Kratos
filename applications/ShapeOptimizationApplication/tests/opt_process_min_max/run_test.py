# Import Kratos core and apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
import os, csv

# =======================================================================================================
# Run analysis
# =======================================================================================================

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

with open('response_combination.csv', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')

    # Skipt header lines
    for itr in range(8):
        next(reader)

    for line in reader:
        combination_value_1 = float(line[1].strip())
        combination_value_2 = float(line[3].strip())
        grad_norm_f1 = line[7].strip()
        grad_norm_f2 = float(line[8].strip())

        TestCase().assertAlmostEqual(combination_value_1, 5.26987E+06, 5)
        TestCase().assertAlmostEqual(combination_value_2, 5.26987E+06, 5)
        TestCase().assertEqual(grad_norm_f1, "-")
        TestCase().assertAlmostEqual(grad_norm_f2, 1.43532E+07, 5)

# Cleaning
kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
kratos_utilities.DeleteDirectoryIfExisting(output_directory)
kratos_utilities.DeleteFileIfExisting("response_combination.csv")
kratos_utilities.DeleteFileIfExisting(os.path.basename(original_directory)+".post.lst")
kratos_utilities.DeleteFileIfExisting(optimization_model_part_name+".post.bin")

# =======================================================================================================