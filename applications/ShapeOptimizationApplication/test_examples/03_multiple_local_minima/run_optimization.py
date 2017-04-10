from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *
import optimization_settings as opt_settings
from math import pi, sin

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
            x = node.X
            z = node.Z
            objective = objective + abs(self.TargetCurveOne(x) - z) * abs(self.TargetCurveTwo(x) - z)
        return objective

    # --------------------------------------------------------------------------        
    def ObjectiveGradient( self, currentDesign ):
        """ Returns the gradient of the objective function """
        sensitivity = {}
        for node in currentDesign.Nodes:
            x = node.X
            z = node.Z
            delta_one = z - self.TargetCurveOne(x)
            delta_two = z - self.TargetCurveTwo(x)
            if abs(delta_one) == 0.0 or abs(delta_two) == 0.0:
                sz = 0.0
            else:
                sz = abs(delta_two) * delta_one / abs(delta_one) + abs(delta_one) * delta_two / abs(delta_two)
            sensitivity[node.Id] = [0.0, 0.0, sz]
        return sensitivity  

    # --------------------------------------------------------------------------        
    def TargetCurveOne( self, x ):
        """ Defines target curve 1 as z=TargetCurveOne(x) """
        if x <= 10.0:
            return 0.0
        elif x<= 30.0:
            return sin(2.0 * pi * (x-10.0) / 40.0)
        else:
            return 0.0

    # --------------------------------------------------------------------------        
    def TargetCurveTwo( self, x ):
        """ Defines target curve 2 as z=TargetCurveTwo(x) """
        if x <= 10.0:
            return 0.0
        elif x<= 30.0:
            return abs(sin(2.0 * pi * (x-10.0) / 20.0))
        else:
            return 0.0       

    # --------------------------------------------------------------------------        

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
