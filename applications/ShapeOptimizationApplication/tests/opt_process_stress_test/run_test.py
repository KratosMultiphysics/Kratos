# Import Kratos core and apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
from KratosMultiphysics.json_output_process import JsonOutputProcess
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
import os, csv

# =======================================================================================================
# Helper functions
# =======================================================================================================
def RunJsonOutputProcess(model):
    parameters = KM.Parameters("""
    {
        "output_variables" : ["SHAPE_CHANGE_X","SHAPE_CHANGE_Y","SHAPE_CHANGE_Z"],
        "output_file_name" : "shape_change_results.json",
        "model_part_name"  : "structure.Parts_structure",
        "time_frequency"   : 0.0
    }
    """)
    process = JsonOutputProcess(model, parameters)
    process.ExecuteInitialize()
    process.ExecuteBeforeSolutionLoop()
    process.ExecuteFinalizeSolutionStep()

def RunJsonCheckProcess(model):
    parameters = KM.Parameters("""
    {
        "check_variables"  : ["SHAPE_CHANGE_X","SHAPE_CHANGE_Y","SHAPE_CHANGE_Z"],
        "input_file_name"  : "shape_change_results.json",
        "model_part_name"  : "structure.Parts_structure",
        "time_frequency"   : 0.0
    }
    """)
    process = FromJsonCheckResultProcess(model, parameters)
    process.ExecuteInitialize()
    process.ExecuteBeforeSolutionLoop()
    process.ExecuteFinalizeSolutionStep()

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
optimization_log_filename = parameters["optimization_settings"]["output"]["optimization_log_filename"].GetString() + ".csv"

# Testing by
# 1) using the "json_output_process" & "json_check_process"
# 2) additionally checking some process output

# Run json check process
# Note that due the adjoint sensitivity process in the structural analysis (CSM) the json_ouput_process and json_check_process cannot be included into the CSM processes
# RunJsonOutputProcess(model)
RunJsonCheckProcess(model)

# Run additional checks
os.chdir(output_directory)

with open(optimization_log_filename, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    last_line = None
    for line in reader:
        if not line:
            continue
        else:
            last_line = line

    resulting_abs_improvement = float(last_line[2].strip())
    TestCase().assertAlmostEqual(resulting_abs_improvement, -62.0931, 4)

    norm_df = float(last_line[4].strip())
    TestCase().assertAlmostEqual(norm_df, 1.93278E+06, 4)

os.chdir(original_directory)

# Cleaning
kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
kratos_utilities.DeleteDirectoryIfExisting(output_directory)
kratos_utilities.DeleteFileIfExisting(os.path.basename(original_directory)+".post.lst")
kratos_utilities.DeleteFileIfExisting(optimization_model_part_name+".time")
kratos_utilities.DeleteFileIfExisting(optimization_model_part_name+".post.bin")

# =======================================================================================================