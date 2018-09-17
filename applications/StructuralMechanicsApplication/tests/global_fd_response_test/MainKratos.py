from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import structural_response_function_factory

with open("strain_energy_response_parameters.json",'r') as parameter_file:
    ProjectParameters = Parameters( parameter_file.read())

model = Model()
response = structural_response_function_factory.CreateResponseFunction("Global_FD", ProjectParameters["kratos_response_settings"], model)
response.RunCalculation(calculate_gradient=True)
