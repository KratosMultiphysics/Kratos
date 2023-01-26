# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
from KratosMultiphysics import Logger

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")

    default_settings = KM.Parameters(
        """{
            "model_part_name"      : "set_model_part_name",
            "cartesian_local_axis" : [[1.0,0.0,0.0],[0.0,1.0,0.0]],
            "update_at_each_step"  : false
        }""");
    process_settings = settings["Parameters"]
    process_settings.ValidateAndAssignDefaults(default_settings)
    computing_model_part = Model[process_settings["model_part_name"].GetString()]

    Logger.PrintInfo("SetCartesianLocalAxesProcess:: ","Setting the oriented local axes...")
    process_settings.RemoveValue("model_part_name")
    return SMA.SetCartesianLocalAxesProcess(computing_model_part, process_settings)




