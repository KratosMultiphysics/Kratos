# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
import csv, os

# Read parameters
with open("parameters.json",'r') as parameter_file:
    parameters = Parameters(parameter_file.read())

# Defining the model_part
optimization_model_part = ModelPart(parameters["optimization_settings"]["design_variables"]["optimization_model_part_name"].GetString())
optimization_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, parameters["optimization_settings"]["design_variables"]["domain_size"].GetInt())

# =======================================================================================================
# Define external analyzer
# =======================================================================================================

# The external analyzer provides a response to constrain the distance of a specific node to a given target
from analyzer_base import AnalyzerBaseClass
class CustomAnalyzer(AnalyzerBaseClass):
    # --------------------------------------------------------------------------------------------------
    def __init__( self ):
        self.constrained_node_id =975
        self.target_x = 1.15655
        self.target_y = 9.93289
        self.target_z = 5.28392

    # --------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):
        if communicator.isRequestingValueOf("distance"):
            communicator.reportValue("distance", self.__CalculateValue())

        if communicator.isRequestingGradientOf("distance"):
            communicator.reportGradient("distance", self.__CalculateGradient())

    # --------------------------------------------------------------------------
    def __CalculateValue( self ):
        constrained_node = optimization_model_part.GetNodes()[self.constrained_node_id]

        distance = [0,0,0]
        distance[0] = constrained_node.X0 - self.target_x
        distance[1] = constrained_node.Y0 - self.target_y
        distance[2] = constrained_node.Z0 - self.target_z

        return distance[0]**2 + distance[1]**2 + distance[2]**2

    # --------------------------------------------------------------------------
    def __CalculateGradient( self ):
        constrained_node = optimization_model_part.GetNodes()[self.constrained_node_id]

        response_gradient = {}
        for node in optimization_model_part.Nodes:

            local_gradient = [0,0,0]

            if node.Id == self.constrained_node_id:
                local_gradient[0] = 2*(constrained_node.X0 - self.target_x)
                local_gradient[1] = 2*(constrained_node.Y0 - self.target_y)
                local_gradient[2] = 2*(constrained_node.Z0 - self.target_z)
                print(local_gradient)
            else:
                local_gradient[0] = 0.0
                local_gradient[1] = 0.0
                local_gradient[2] = 0.0

            response_gradient[node.Id] = local_gradient

        return response_gradient

# =======================================================================================================
# Perform optimization
# =======================================================================================================

# Create optimizer and perform optimization
import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], optimization_model_part, CustomAnalyzer())
optimizer.Optimize()

# =======================================================================================================
# Test results and clean directory
# =======================================================================================================
output_directory = parameters["optimization_settings"]["output"]["output_directory"].GetString()
response_log_filename = parameters["optimization_settings"]["output"]["response_log_filename"].GetString() + ".csv"
optimization_model_part_name = parameters["optimization_settings"]["design_variables"]["optimization_model_part_name"].GetString()

# Testing
original_directory = os.getcwd()
os.chdir(output_directory)

with open(response_log_filename, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    last_line = None
    for line in reader:
        if not line:
            continue
        else:
            last_line = line

    resulting_improvement = float(last_line[2].strip())
    resulting_constraint_value = float(last_line[4].strip())

    # # Check against specifications
    TestCase().assertAlmostEqual(resulting_improvement, -19.201365, 2)
    TestCase().assertAlmostEqual(resulting_constraint_value, 0.015641, 4)

os.chdir(original_directory)

# Cleaning
kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
kratos_utilities.DeleteDirectoryIfExisting(output_directory)
kratos_utilities.DeleteFileIfExisting(os.path.basename(original_directory)+".post.lst")
kratos_utilities.DeleteFileIfExisting(optimization_model_part_name+".time")

# =======================================================================================================