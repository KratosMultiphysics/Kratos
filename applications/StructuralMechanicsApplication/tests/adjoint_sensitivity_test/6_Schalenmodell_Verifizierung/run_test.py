from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import structural_mechanics_analysis

# Create the primal solver
with open("rectangular_plate_parameters.json",'r') as parameter_file:
    ProjectParametersPrimal = Parameters( parameter_file.read())
primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersPrimal)

# Create the adjoint solver
with open("rectangular_plate_adjoint_parameters.json",'r') as parameter_file:
    ProjectParametersAdjoint = Parameters( parameter_file.read())
adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)

primal_analysis.Run()
adjoint_analysis.Run()