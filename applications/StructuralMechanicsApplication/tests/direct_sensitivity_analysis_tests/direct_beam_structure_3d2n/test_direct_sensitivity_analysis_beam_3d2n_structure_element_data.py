from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
#import KratosMultiphysics
#from KratosMultiphysics.StructuralMechanicsApplication import *
#import structural_mechanics_analysis
#import KratosMultiphysics.kratos_utilities as kratos_utilities
#import KratosMultiphysics.ExternalSolversApplication
#import direct_sensitivity_analysis 

from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_mechanics_analysis
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.ExternalSolversApplication
import direct_sensitivity_analysis

try:
    from KratosMultiphysics.HDF5Application import *
    has_hdf5_application = True
except ImportError:
    has_hdf5_application = False


with open("direct_sensitivity_beam_parameters_element_data.json",'r') as parameter_file:
    ProjectParameters = KratosMultiphysics.Parameters( parameter_file.read())

model = Model()

direct_sensitivity_analysis = direct_sensitivity_analysis.DirectSensitivityAnalysis(ProjectParameters,model)
calculate_gradient = True
direct_sensitivity_analysis.RunCalculation(calculate_gradient)
print("Simulation finished")
