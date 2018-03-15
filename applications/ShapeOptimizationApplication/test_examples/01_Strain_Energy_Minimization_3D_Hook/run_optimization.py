from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# For time measures
import time as timer

# Read parameters
with open("ProjectParameters.json",'r') as parameter_file:
    project_parameters = Parameters(parameter_file.read())

# Defining the model_part
optimization_model_part = ModelPart(project_parameters["optimization_settings"]["design_variables"]["optimization_model_part_name"].GetString())
optimization_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, project_parameters["optimization_settings"]["design_variables"]["domain_size"].GetInt())

# Create optimizer and perform optimization
import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer(project_parameters, optimization_model_part)
optimizer.optimize()