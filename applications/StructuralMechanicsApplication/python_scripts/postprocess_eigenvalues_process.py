# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as KSM

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")

    process_settings = settings["Parameters"]
    __BackwardsCompatibilityHelper(process_settings)
    return KSM.PostprocessEigenvaluesProcess(Model, process_settings)

def __BackwardsCompatibilityHelper(process_settings):
    # Check if "computing_model_part_name" is provided
    if process_settings.Has("computing_model_part_name"):
        KratosMultiphysics.Logger.PrintWarning("'computing_model_part_name' is deprecated. Use 'model_part_name' instead.")
        process_settings.AddEmptyValue("model_part_name").SetString(process_settings["computing_model_part_name"].GetString())
        process_settings.RemoveValue("computing_model_part_name")
    # Remove the old "help" field
    process_settings.RemoveValue("help")

