from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as KSM

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")

    if settings.Has("computing_model_part_name"):
        computing_model_part = Model[settings["computing_model_part_name"].GetString()]
    else: # using default name
        computing_model_part = Model["Structure.computing_domain"]

    settings.RemoveValue("help")

    return KSM.PostprocessEigenvaluesProcess(computing_model_part, settings)
