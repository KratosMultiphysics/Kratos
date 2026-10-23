# Import Kratos core and apps
import KratosMultiphysics as KM

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.ShapeOptimizationApplication.analyzers.analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess

import os, csv

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

model = KM.Model()

class CustomAnalyzer(AnalyzerBaseClass):
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):
        if communicator.isRequestingValueOf("reciprocal_r_squared_sum"):
            value = 0.0
            for node in current_design.Nodes:
                r_squared = node.X**2 + node.Y**2
                value += r_squared
            value = 1/value
            communicator.reportValue("reciprocal_r_squared_sum", value)

        if communicator.isRequestingGradientOf("reciprocal_r_squared_sum"):
            gradient = {}
            for node in current_design.Nodes:
                r_squared = node.X**2 + node.Y**2
                gradient[node.Id] = [-2*node.X/r_squared**2, -2*node.Y/r_squared**2, 0.0]
            communicator.reportGradient("reciprocal_r_squared_sum", gradient)

# Create optimizer and perform optimization
optimizer = optimizer_factory.Create(model, parameters["optimization_settings"], CustomAnalyzer())
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
    TestCase().assertEqual(resulting_iteration, 1)

# # write json output
# output_process = KM.JsonOutputProcess(model, KM.Parameters(
#     """{
#         "output_variables" : ["SHAPE_UPDATE_X","SHAPE_UPDATE_Y","SHAPE_UPDATE_Z"],
#         "output_file_name" : "shape_update_results.json",
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
        "check_variables"  : ["SHAPE_UPDATE_X","SHAPE_UPDATE_Y","SHAPE_UPDATE_Z"],
        "input_file_name"  : "shape_update_results.json",
        "model_part_name"  : \""""+optimization_model_part_name+"""\",
        "time_frequency"   : 0.0
    }"""))
check_process.ExecuteInitialize()
check_process.ExecuteBeforeSolutionLoop()
check_process.ExecuteInitializeSolutionStep()
check_process.ExecuteFinalizeSolutionStep()
check_process.ExecuteFinalize()

# =======================================================================================================