# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters = Parameters(parameter_file.read())

# Defining the model_part
optimization_model_part = ModelPart(parameters["optimization_settings"]["design_variables"]["optimization_model_part_name"].GetString())
optimization_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, parameters["optimization_settings"]["design_variables"]["domain_size"].GetInt())

# Create optimizer and perform optimization
import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], optimization_model_part)
optimizer.Optimize()