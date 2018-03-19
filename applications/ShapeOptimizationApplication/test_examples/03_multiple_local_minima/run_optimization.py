# This test example is from M. Hojjat, E. Stavropoulou,
# K.-U. Bletzinger, The Vertex Morphing method for node-based
# shape optimization, Comput. Methods Appl. Mech. Engrg. 268
# (2014) 494-513.
#
# The target curves are illustrated below.
#
#                    z=1
#                **  +++  **    target curve 1 (+)
#              *  + *   * +  *  target curve 2 (*)
#             * +    * *    + *
#  |--> x     +               +     z=0
#  __________*________*________*__________
#  |-- 10 ---|-- 10 --|-- 10 --|--- 10 --|
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
            x = node.X
            z = node.Z
            objective = objective + abs(self.__TargetCurveOne(x) - z) * abs(self.__TargetCurveTwo(x) - z)
        return objective

    # --------------------------------------------------------------------------------------------------
    def __ObjectiveGradient(self, current_design):
        """ Returns the gradient of the objective function """
        sensitivity = {}
        for node in current_design.Nodes:
            x = node.X
            z = node.Z
            delta_one = z - self.__TargetCurveOne(x)
            delta_two = z - self.__TargetCurveTwo(x)
            if abs(delta_one) == 0.0 or abs(delta_two) == 0.0:
                sz = 0.0
            else:
                sz = abs(delta_two) * delta_one / abs(delta_one) + abs(delta_one) * delta_two / abs(delta_two)
            sensitivity[node.Id] = [0.0, 0.0, sz]
        return sensitivity

    # --------------------------------------------------------------------------------------------------
    def __TargetCurveOne(self, x):
        """ Defines target curve 1 as z=__TargetCurveOne(x) """
        if x <= 10.0:
            return 0.0
        elif x<= 30.0:
            return sin(2.0 * pi * (x-10.0) / 40.0)
        else:
            return 0.0

    # --------------------------------------------------------------------------------------------------
    def __TargetCurveTwo(self, x):
        # Defines target curve 2 as z=TargetCurveTwo(x)
        if x <= 10.0:
            return 0.0
        elif x<= 30.0:
            return abs(sin(2.0 * pi * (x-10.0) / 20.0))
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