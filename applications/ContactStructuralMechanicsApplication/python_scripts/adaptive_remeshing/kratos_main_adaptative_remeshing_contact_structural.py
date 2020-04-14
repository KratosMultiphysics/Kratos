from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as KM

from KratosMultiphysics.ContactStructuralMechanicsApplication.adaptive_remeshing.adaptative_remeshing_contact_structural_mechanics_analysis import AdaptativeRemeshingContactStructuralMechanicsAnalysis

## Import define_output
with open("ProjectParameters.json",'r') as parameter_file:
    ProjectParameters = KM.Parameters(parameter_file.read())

# Creating the test
model = KM.Model()
analysis = AdaptativeRemeshingContactStructuralMechanicsAnalysis(model, ProjectParameters)
analysis.Run()
