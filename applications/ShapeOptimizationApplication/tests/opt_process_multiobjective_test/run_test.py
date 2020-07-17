# Import Kratos core and apps
import KratosMultiphysics as KM

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.ShapeOptimizationApplication.analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
import os, csv

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

model = KM.Model()

# =======================================================================================================
# Define external analyzer
# =======================================================================================================
class CustomAnalyzer(AnalyzerBaseClass):
    # --------------------------------------------------------------------------------------------------
    def __init__( self ):
        self.measurement_node_id = 5

    # --------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):
        if communicator.isRequestingValueOf("length"):
            communicator.reportValue("length", self.__CalculateValue(current_design))

        if communicator.isRequestingGradientOf("length"):
            communicator.reportGradient("length", self.__CalculateGradient(current_design))

    # --------------------------------------------------------------------------
    def __CalculateValue( self, current_design ):
        node = current_design.GetNodes()[self.measurement_node_id]
        return node.X

    # --------------------------------------------------------------------------
    def __CalculateGradient( self, current_design ):
        node = current_design.GetNodes()[self.measurement_node_id]

        response_gradient = {}
        for node in current_design.Nodes:
            local_gradient = [0,0,0]

            if node.Id == self.measurement_node_id:
                local_gradient[0] = 1.0
                local_gradient[1] = 0.0
                local_gradient[2] = 0.0

            response_gradient[node.Id] = local_gradient

        return response_gradient

# =======================================================================================================
# Perform optimization
# =======================================================================================================

# Create optimizer and perform optimization
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model, CustomAnalyzer())
optimizer.Optimize()

# =======================================================================================================
# Test results and clean directory
# =======================================================================================================
# Testing by
# 1) using the "json_output_process" & "json_check_process" within the structural analysis

# 2) checking values of response function combinations
response_combination_filename = "response_combination.csv"

with open(response_combination_filename, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    last_line = None
    for line in reader:
        if not line:
            continue
        else:
            last_line = line

    f1_value = float(last_line[2].strip())
    f3_value = float(last_line[4].strip())
    f4_value = float(last_line[5].strip())

    TestCase().assertAlmostEqual(f1_value, 5.03120E-01, 4)
    TestCase().assertAlmostEqual(f3_value, 2.17232E+02, 4)
    TestCase().assertAlmostEqual(f4_value, 1.20762E+00, 4)

# 3) additionally checking some optimization output
original_directory = os.getcwd()
output_directory = parameters["optimization_settings"]["output"]["output_directory"].GetString()
optimization_model_part_name = parameters["optimization_settings"]["model_settings"]["model_part_name"].GetString()
optimization_log_filename = parameters["optimization_settings"]["output"]["optimization_log_filename"].GetString() + ".csv"

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
    TestCase().assertAlmostEqual(resulting_abs_improvement, -1.72389E+01, 4)

os.chdir(original_directory)

# =======================================================================================================