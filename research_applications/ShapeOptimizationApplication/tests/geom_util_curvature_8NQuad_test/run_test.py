# Import Kratos core and apps
import KratosMultiphysics as Kratos

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.ShapeOptimizationApplication.analyzers.analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.KratosUnittest import TestCase
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess
# from KratosMultiphysics.json_output_process import JsonOutputProcess

import os, csv

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters = Kratos.Parameters(parameter_file.read())

model = Kratos.Model()

class CustomAnalyzer(AnalyzerBaseClass):
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):
        if communicator.isRequestingValueOf("no_objective"):
            value = 0.0
            communicator.reportValue("no_objective", value)

        if communicator.isRequestingGradientOf("no_objective"):
            gradient = {}
            for node in current_design.Nodes:
                gradient[node.Id] = [0.0, 0.0, 0.0]
            communicator.reportGradient("no_objective", gradient)


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
# output_process = JsonOutputProcess(model, Kratos.Parameters(
#     """{
#         "output_variables" : ["GAUSSIAN_CURVATURE"],
#         "output_file_name" : "gaussian_curvature_results.json",
#         "model_part_name"  : \""""+optimization_model_part_name+"""\",
#         "time_frequency"   : 0.0
#     }"""))

# output_process.ExecuteInitialize()
# output_process.ExecuteBeforeSolutionLoop()
# output_process.ExecuteInitializeSolutionStep()
# output_process.ExecuteFinalizeSolutionStep()
# output_process.ExecuteFinalize()

check_process = FromJsonCheckResultProcess(model, Kratos.Parameters(
    """{
        "check_variables"  : ["GAUSSIAN_CURVATURE"],
        "input_file_name"  : "gaussian_curvature_results.json",
        "model_part_name"  : \""""+optimization_model_part_name+"""\",
        "time_frequency"   : 0.0
    }"""))
check_process.ExecuteInitialize()
check_process.ExecuteBeforeSolutionLoop()
check_process.ExecuteInitializeSolutionStep()
check_process.ExecuteFinalizeSolutionStep()
check_process.ExecuteFinalize()

# =======================================================================================================