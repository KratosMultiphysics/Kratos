from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *
import optimization_settings as opt_settings
import sys

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

import optimization_settings as optimizationSettings
import optimizer_factory as optimizerFactory

# ======================================================================================================================================
# Solver preparation
# ======================================================================================================================================

class externalAnalyzer( optimizerFactory.analyzerBaseClass ):
    
    # --------------------------------------------------------------------------
    def analyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):
        if communicator.isRequestingFunctionValueOf("targetDeviation"): 
            communicator.reportFunctionValue("targetDeviation", self.ObjectiveFunction(currentDesign))    
        if communicator.isRequestingGradientOf("targetDeviation"): 
            communicator.reportGradient("targetDeviation", self.ObjectiveGradient(currentDesign))   

    # --------------------------------------------------------------------------
    def ObjectiveFunction( self, currentDesign ):
        """ Returns the objective function to be minimized """
        objective = 0.0
        for node in currentDesign.Nodes:
            objective = objective + abs(self.TentFunction(node.X) - node.Z)
        return objective

    # --------------------------------------------------------------------------
    def ObjectiveGradient( self, currentDesign ):
        """ Returns the gradient of the objective function """
        sensitivity = dict()
        for node in currentDesign.Nodes:
            delta = node.Z - self.TentFunction(node.X)
            if abs(delta) == 0.0:
                sz = 0.0
            else:
                sz = delta / abs(delta)
            sensitivity[node.Id] = [0.0, 0.0, sz]
        return sensitivity

    # --------------------------------------------------------------------------
    def TentFunction( self, x ):
        """ Defines the target curve z=TentFunction(x) """
        if x <= 15.0:
            return 0.0
        elif x<= 20.0:
            return (x - 15.0) / 5.0
        elif x <= 25.0:
            return 1.0 - (x - 20.0) / 5.0
        else:
            return 0.0


# ======================================================================================================================================
# Optimization part
# ======================================================================================================================================

inputModelPart = ModelPart(optimizationSettings.input_model_part_name)
newAnalyzer = externalAnalyzer()

optimizer = optimizerFactory.CreateOptimizer( inputModelPart, optimizationSettings )
optimizer.importAnalyzer( newAnalyzer )
optimizer.importModelPart()

optimizer.optimize()

# ======================================================================================================================================
