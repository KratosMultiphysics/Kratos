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

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
import csv, os

# =======================================================================================================
# Define external analyzer
# =======================================================================================================

from analyzer_base import AnalyzerBaseClass
class CustomAnalyzer(AnalyzerBaseClass):

    # --------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):
        if communicator.isRequestingValueOf("targetDeviation"):
            communicator.reportValue("targetDeviation", self.__ObjectiveFunction(current_design))

        if communicator.isRequestingGradientOf("targetDeviation"):
            communicator.reportGradient("targetDeviation", self.__ObjectiveGradient(current_design))

    # --------------------------------------------------------------------------------------------------
    def __ObjectiveFunction(self, current_design):
        """ Returns the objective function to be minimized """
        objective = 0.0
        for node in current_design.Nodes:
            objective = objective + abs(self.__TentFunction(node.X) - node.Z)
        return objective

    # --------------------------------------------------------------------------------------------------
    def __ObjectiveGradient(self, current_design):
        """ Returns the gradient of the objective function """
        sensitivity = dict()
        for node in current_design.Nodes:
            delta = node.Z - self.__TentFunction(node.X)
            if abs(delta) == 0.0:
                sz = 0.0
            else:
                sz = delta / abs(delta)
            sensitivity[node.Id] = [0.0, 0.0, sz]
        return sensitivity

    # --------------------------------------------------------------------------------------------------
    @staticmethod
    def __TentFunction(x):
        """ Defines the target curve z=__TentFunction(x) """
        if x <= 15.0:
            return 0.0
        elif x<= 20.0:
            return (x - 15.0) / 5.0
        elif x <= 25.0:
            return 1.0 - (x - 20.0) / 5.0
        else:
            return 0.0

# =======================================================================================================
# Perform optimization
# =======================================================================================================

with open("parameters.json",'r') as parameter_file:
    parameters = Parameters(parameter_file.read())

optimization_model_part = ModelPart(parameters["optimization_settings"]["design_variables"]["optimization_model_part_name"].GetString())
optimization_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, parameters["optimization_settings"]["design_variables"]["domain_size"].GetInt())

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

    resulting_optimization_iterations = int(last_line[0].strip())
    resulting_improvement = float(last_line[2].strip())

    # Check against specifications
    TestCase().assertEqual(resulting_optimization_iterations, 124)
    TestCase().assertAlmostEqual(resulting_improvement, -90.228425, 5)

os.chdir(original_directory)

# Cleaning
kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
kratos_utilities.DeleteDirectoryIfExisting(output_directory)
kratos_utilities.DeleteFileIfExisting(os.path.basename(original_directory)+".post.lst")
kratos_utilities.DeleteFileIfExisting(optimization_model_part_name+".time")

# =======================================================================================================
