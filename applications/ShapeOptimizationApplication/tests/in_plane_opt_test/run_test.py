# Import Kratos core and apps
import KratosMultiphysics as KM

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.ShapeOptimizationApplication.analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess
from KratosMultiphysics.json_output_process import JsonOutputProcess

import os, csv

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

model = KM.Model()

class CustomAnalyzer(AnalyzerBaseClass):
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):
        if communicator.isRequestingValueOf("x_squared_sum"):
            value = 0.0
            for node in current_design.Nodes:
                value += node.X**2
            communicator.reportValue("x_squared_sum", value)

        if communicator.isRequestingGradientOf("x_squared_sum"):
            gradient = {}
            for node in current_design.Nodes:
                gradient[node.Id] = [2*node.X, 0.0, 0.0]
            communicator.reportGradient("x_squared_sum", gradient)

# Create optimizer and perform optimization
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model, CustomAnalyzer())
optimizer.Optimize()

# =======================================================================================================
# Test results and clean directory
# =======================================================================================================
output_directory = parameters["optimization_settings"]["output"]["output_directory"].GetString()
optimization_model_part_name = parameters["optimization_settings"]["model_settings"]["model_part_name"].GetString()
optimization_log_filename = parameters["optimization_settings"]["output"]["optimization_log_filename"].GetString() + ".csv"

# Testing by
# 1) some values from csv output
# 2) using the "json_output_process" & "json_check_process"

with open(os.path.join(output_directory, optimization_log_filename), 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    last_line = None
    for line in reader:
        if not line:
            continue
        else:
            last_line = line

    resulting_iteration = float(last_line[0].strip())
    resulting_abs_improvement = float(last_line[2].strip())
    resulting_constraint = float(last_line[4].strip())

    # Check against specifications
    TestCase().assertEqual(resulting_iteration, 5)
    TestCase().assertAlmostEqual(resulting_abs_improvement, -3.41106E+01, 3)

# # write json output
# output_process = JsonOutputProcess(model, KM.Parameters(
#     """{
#         "output_variables" : ["SHAPE_CHANGE_X","SHAPE_CHANGE_Y","SHAPE_CHANGE_Z"],
#         "output_file_name" : "shape_change_results.json",
#         "model_part_name"  : \""""+optimization_model_part_name+"""\",
#         "time_frequency"   : 0.0
#     }"""))

# output_process.ExecuteInitialize()
# output_process.ExecuteBeforeSolutionLoop()
# output_process.ExecuteInitializeSolutionStep()
# output_process.ExecuteFinalizeSolutionStep()
# output_process.ExecuteFinalize()

check_process = FromJsonCheckResultProcess(model, KM.Parameters(
    """{
        "check_variables"  : ["SHAPE_CHANGE_X","SHAPE_CHANGE_Y","SHAPE_CHANGE_Z"],
        "input_file_name"  : "shape_change_results.json",
        "model_part_name"  : \""""+optimization_model_part_name+"""\",
        "time_frequency"   : 0.0
    }"""))
check_process.ExecuteInitialize()
check_process.ExecuteBeforeSolutionLoop()
check_process.ExecuteInitializeSolutionStep()
check_process.ExecuteFinalizeSolutionStep()
check_process.ExecuteFinalize()

# Cleaning
kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
kratos_utilities.DeleteDirectoryIfExisting(output_directory)
kratos_utilities.DeleteFileIfExisting(os.path.basename(os.getcwd())+".post.lst")

# =======================================================================================================