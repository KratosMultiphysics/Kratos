from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as KSM

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")

    process_settings = settings["Parameters"]

    if process_settings.Has("model_part_name"):
        computing_model_part = model[process_settings["model_part_name"].GetString()]
    else: # using default name
        computing_model_part = model["Structure"]

    process_settings.RemoveValue("model_part_name")
    process_settings.RemoveValue("help")

    return KSM.GeneralizedInfluenceFunctionsProcess(computing_model_part, process_settings)










