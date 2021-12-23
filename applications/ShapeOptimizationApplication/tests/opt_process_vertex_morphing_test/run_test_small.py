# This test example is from M. Hojjat, E. Stavropoulou,
# K.-U. Bletzinger, The Vertex Morphing method for node-based
# shape optimization, Comput. Methods Appl. Mech. Engrg. 268
# (2014) 494-513.
#
# The target curve is the tent function illustrated below.
#
#                    z=1
#                     /\
#                    /  \
#                   /    \
#  |--> x          /      \           z=0
#  _______________/        \_______________
#  |----- 15 -----|-- 10 --|----- 15 -----|
#
#
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
    from .tent_analyzer import CustomAnalyzer
except ImportError:
    from tent_analyzer import CustomAnalyzer

# =======================================================================================================
# Perform optimization
# =======================================================================================================

with open("parameters.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

# reduce number of iterations to conform to the maximum execution time in CI
parameters["optimization_settings"]["optimization_algorithm"]["max_iterations"].SetInt(5)
parameters["optimization_settings"]["optimization_algorithm"]["line_search"]["step_size"].SetDouble(0.1)

# activate design output for testing
parameters["optimization_settings"]["output"]["design_output_mode"].SetString("write_design_surface")

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
    TestCase().assertEqual(resulting_optimization_iterations, 5)
    TestCase().assertAlmostEqual(resulting_improvement, -26.2433, 4)

# Testing of design output
output_file_name = os.path.join(output_directory, "design_surface_0_5.vtk")
reference_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ref_design_surface_0_5.vtk")

result_check_settings = KM.Parameters("""{
        "reference_file_name"   : \""""+reference_file_name.replace("\\", "\\\\")+"""\",
        "output_file_name"      : \""""+output_file_name.replace("\\", "\\\\")+"""\",
        "remove_output_file"    : false,
        "comparison_type"       : "vtk",
        "tolerance"             : 1e-6,
        "relative_tolerance"    : 1e-9,
        "dimension"             : 3
    }""")

results_check_process = CompareTwoFilesCheckProcess(result_check_settings)
results_check_process.Execute()

# =======================================================================================================
