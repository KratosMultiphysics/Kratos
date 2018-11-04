# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
import KratosMultiphysics.kratos_utilities as kratos_utilities
import os

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

# Testing is done using the "json_output_process" & "json_check_process" within the structural analysis

# Cleaning
original_directory = os.getcwd()
output_directory = parameters["optimization_settings"]["output"]["output_directory"].GetString()
optimization_model_part_name = parameters["optimization_settings"]["model_settings"]["model_part_name"].GetString()
kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
kratos_utilities.DeleteDirectoryIfExisting(output_directory)
kratos_utilities.DeleteFileIfExisting(os.path.basename(original_directory)+".post.lst")
kratos_utilities.DeleteFileIfExisting(optimization_model_part_name+".time")

# =======================================================================================================