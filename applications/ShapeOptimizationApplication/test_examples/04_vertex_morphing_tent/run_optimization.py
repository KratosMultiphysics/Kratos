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
from math import pi, sin

# =======================================================================================================
# Define external analyzer
# =======================================================================================================

from analyzer_base import AnalyzerBaseClass
class CustomAnalyzer( AnalyzerBaseClass ):

    # --------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator( self, currentDesign, OptimizationIteration, Communicator ):
        if Communicator.isRequestingValueOf("targetDeviation"):
            Communicator.reportValue("targetDeviation", self.__ObjectiveFunction(currentDesign))

        if Communicator.isRequestingGradientOf("targetDeviation"):
            Communicator.reportGradient("targetDeviation", self.__ObjectiveGradient(currentDesign))

    # --------------------------------------------------------------------------------------------------
    def __ObjectiveFunction( self, currentDesign ):
        """ Returns the objective function to be minimized """
        objective = 0.0
        for node in currentDesign.Nodes:
            objective = objective + abs(self.__TentFunction(node.X) - node.Z)
        return objective

    # --------------------------------------------------------------------------------------------------
    def __ObjectiveGradient( self, currentDesign ):
        """ Returns the gradient of the objective function """
        sensitivity = dict()
        for node in currentDesign.Nodes:
            delta = node.Z - self.__TentFunction(node.X)
            if abs(delta) == 0.0:
                sz = 0.0
            else:
                sz = delta / abs(delta)
            sensitivity[node.Id] = [0.0, 0.0, sz]
        return sensitivity

    # --------------------------------------------------------------------------------------------------
    def __TentFunction( self, x ):
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

with open("ProjectParameters.json",'r') as parameter_file:
    parameters = Parameters(parameter_file.read())

optimization_model_part = ModelPart(parameters["optimization_settings"]["design_variables"]["optimization_model_part_name"].GetString())
optimization_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, parameters["optimization_settings"]["design_variables"]["domain_size"].GetInt())

import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer(parameters, optimization_model_part, CustomAnalyzer())
optimizer.Optimize()

# =======================================================================================================