# Import Kratos core and apps
import KratosMultiphysics as KM

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess
from KratosMultiphysics.json_output_process import JsonOutputProcess

import os, csv

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

model = KM.Model()

# Create optimizer and perform optimization
optimizer = optimizer_factory.Create(model, parameters["optimization_settings"])
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

# TODO: finite difference sensitivity testing
#       - implement area derivatives
#       - fd sensitivities can only be computed for the infeasible nodes ("volume nodes") which are forming a pond
#       - feasible nodes don't have any sensitivity

objective_reference_result = [2.11349E+01, 5.28372E+00, 0.00000E+00]

with open(os.path.join(output_directory, optimization_log_filename), 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    last_line = None
    resulting_objective = []
    i = 0
    for line in reader:
        if not line:
            continue
        else:
            if i != 0:
                resulting_objective.append(float(line[1].strip()))
            last_line = line
            i += 1

    resulting_iteration = float(last_line[0].strip())
    resulting_abs_improvement = float(last_line[2].strip())

    # Check against specifications
    TestCase().assertEqual(resulting_iteration, 3)
    TestCase().assertAlmostEqual(resulting_objective, objective_reference_result, 3)

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

# =======================================================================================================