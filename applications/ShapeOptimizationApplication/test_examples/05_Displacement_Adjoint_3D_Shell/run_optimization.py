from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *
import optimizer_factory # from ShapeOptimizationApplication

# For time measures
import time as timer

# ======================================================================================================================================
# Model part & solver
# ======================================================================================================================================

with open("optimization_parameters.json",'r') as parameter_file:
    parameters = Parameters( parameter_file.read())

# Defining the model_part
optimization_model_part = ModelPart(parameters["optimization_settings"]["design_variables"]["optimization_model_part_name"].GetString())
optimization_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, parameters["optimization_settings"]["design_variables"]["domain_size"].GetInt())

# Create optimizer and perform optimization
optimizer = optimizer_factory.CreateOptimizer(parameters, optimization_model_part)
optimizer.Optimize()