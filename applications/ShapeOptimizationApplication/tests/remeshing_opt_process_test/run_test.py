# Import Kratos core and apps
import KratosMultiphysics as KM
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication.analyzers.analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
import csv, os

try:
    from KratosMultiphysics.ShapeOptimizationApplication.analyzers.custom_analyzer_y import CustomAnalyzer
except ImportError:
    from .custom_analyzer_y import CustomAnalyzer

# =======================================================================================================
# Perform optimization
# =======================================================================================================

with open("parameters.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

model = KM.Model()

optimizer = optimizer_factory.Create(model, parameters["optimization_settings"], CustomAnalyzer())
optimizer.Optimize()

# =======================================================================================================
# Test results and clean directory
# =======================================================================================================
output_directory = parameters["optimization_settings"]["output"]["output_directory"].GetString()
optimization_log_filename = parameters["optimization_settings"]["output"]["optimization_log_filename"].GetString() + ".csv"

# Testing of optimization progress
with open(os.path.join(output_directory, optimization_log_filename), 'r') as csvfile:
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
    TestCase().assertEqual(resulting_optimization_iterations, 2)
    TestCase().assertAlmostEqual(resulting_improvement, 1.42057E+01, 4)

# check if remeshing has happened
mp = model.GetModelPart(parameters["optimization_settings"]["model_settings"]["model_part_name"].GetString())

TestCase().assertNotEqual(mp.NumberOfNodes(), 422) # initial mesh had 422 nodes
TestCase().assertEqual(mp.NumberOfNodes(), 525) # after remeshing in step 1 it has 525 nodes

# =======================================================================================================
