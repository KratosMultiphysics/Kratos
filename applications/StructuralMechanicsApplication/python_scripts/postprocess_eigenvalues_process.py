from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as KSM

# Other imports
import os

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")

    process_settings = settings["Parameters"]

    folder_settings = KratosMultiphysics.Parameters("""{
        "folder_name"                 : "EigenResults",
        "save_output_files_in_folder" : true
    }""")

    process_settings.AddMissingParameters(folder_settings)

    if process_settings["save_output_files_in_folder"].GetBool():
        folder_name = process_settings["folder_name"].GetString()
        kratos_utils.DeleteDirectoryIfExisting(folder_name) # make sure to remove old results
        os.mkdir(folder_name)

    if process_settings.Has("computing_model_part_name"):
        computing_model_part = Model[process_settings["computing_model_part_name"].GetString()]
    else: # using default name
        computing_model_part = Model["Structure"]

    process_settings.RemoveValue("computing_model_part_name")
    process_settings.RemoveValue("help")

    return KSM.PostprocessEigenvaluesProcess(computing_model_part, process_settings)
