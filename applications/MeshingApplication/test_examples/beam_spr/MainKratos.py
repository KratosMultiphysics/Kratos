from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.MeshingApplication as MA
import KratosMultiphysics.StructuralMechanicsApplication as SMA

import adaptative_remeshing_structural_mechanics_analysis as Analysis

## Import define_output
with open("ProjectParameters.json",'r') as parameter_file:
    ProjectParameters = KM.Parameters(parameter_file.read())

# Creating the test
model = KM.Model()
analysis = Analysis.AdaptativeRemeshingStructuralMechanicsAnalysis(model, ProjectParameters)
KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.INFO)
analysis.Run()
