from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *
import optimization_settings as settings
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

design_surface = ModelPart("design_surface")

def TentFunction(x):
    """ Defines the target curve z=TentFunction(x) """
    if x <= 15.0:
        return 0.0
    elif x<= 20.0:
        return (x - 15.0) / 5.0
    elif x <= 25.0:
        return 1.0 - (x - 20.0) / 5.0
    else:
        return 0.0

def ObjectiveFunction(X,opt_iter):
    """ Returns the objective function to be minimized """
    objective = 0.0
    for node in design_surface.Nodes:
        x = node.X + X[node.Id][0]
        z = node.Z + X[node.Id][2]
        objective = objective + abs(TentFunction(x) - z)
    return objective

def ObjectiveGradient(X,opt_iter):
    """ Returns the gradient of the objective function """
    sensitivity = dict()
    for node in design_surface.Nodes:
        x = node.X + X[node.Id][0]
        z = node.Z + X[node.Id][2]
        delta = z - TentFunction(x)
        if abs(delta) == 0.0:
            sz = 0.0
        else:
            sz = delta / abs(delta)
        sensitivity[node.Id] = [0.0, 0.0, sz]
    return sensitivity

def Analyzer(X, controls, iterator, response):
    """ Performs the analysis in the optimization loop """
    response["1"]["func"] = ObjectiveFunction(X,iterator)
    response["1"]["grad"] = ObjectiveGradient(X,iterator)

import optimizer_factory

optimizer = optimizer_factory.CreateOptimizer(design_surface,settings.KratosShapeSettings,Analyzer)

optimizer.optimize()
