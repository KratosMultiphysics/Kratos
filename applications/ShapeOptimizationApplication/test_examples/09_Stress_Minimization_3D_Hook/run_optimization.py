# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *
import optimizer_factory # from ShapeOptimizationApplication

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters = Parameters(parameter_file.read())

# Defining the model_part
model = Model()

# Create optimizer and perform optimization
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model)
optimizer.Optimize()