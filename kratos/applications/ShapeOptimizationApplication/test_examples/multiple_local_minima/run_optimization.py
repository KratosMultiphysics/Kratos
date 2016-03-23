from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *
import optimization_settings as settings
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

design_surface = ModelPart("design_surface")

def TargetCurveOne(x):
    """ Defines target curve 1 as z=TargetCurveOne(x) """
    if x <= 10.0:
        return 0.0
    elif x<= 30.0:
        return sin(2.0 * pi * (x-10.0) / 40.0)
    else:
        return 0.0

def TargetCurveTwo(x):
    """ Defines target curve 2 as z=TargetCurveTwo(x) """
    if x <= 10.0:
        return 0.0
    elif x<= 30.0:
        return abs(sin(2.0 * pi * (x-10.0) / 20.0))
    else:
        return 0.0

def ObjectiveFunction(X,opt_iter):
    """ Returns the objective function to be minimized """
    objective = 0.0
    for node in design_surface.Nodes:
        x = node.X
        z = node.Z
        objective = objective + abs(TargetCurveOne(x) - z) \
                    * abs(TargetCurveTwo(x) - z)
    return objective

def ObjectiveGradient(X,opt_iter):
    """ Returns the gradient of the objective function """
    sensitivity = dict()
    for node in design_surface.Nodes:
        x = node.X
        z = node.Z
        delta_one = z - TargetCurveOne(x)
        delta_two = z - TargetCurveTwo(x)
        if abs(delta_one) == 0.0 or abs(delta_two) == 0.0:
            sz = 0.0
        else:
            sz = abs(delta_two) * delta_one / abs(delta_one) \
                 + abs(delta_one) * delta_two / abs(delta_two)
        sensitivity[node.Id] = [0.0, 0.0, sz]
    return sensitivity

def Analyzer(X, controls, iterator, response):
    """ Performs the analysis in the optimization loop """
    if controls["1"]["calc_func"]:
        response["1"]["func"] = ObjectiveFunction(X,iterator)
    if controls["1"]["calc_grad"]:
        response["1"]["grad"] = ObjectiveGradient(X,iterator)

import optimizer_factory

optimizer = optimizer_factory.CreateOptimizer(design_surface,settings.KratosShapeSettings,Analyzer)

optimizer.optimize()
