# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
from KratosMultiphysics import Logger


def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")

    default_settings = KM.Parameters(
        """{
            "model_part_name"               : "set_model_part_name",
            "cylindrical_generatrix_axis"   : [0.0,0.0,1.0],
            "cylindrical_generatrix_point"  : [0.0,0.0,0.0]
        }""")
    process_settings = settings["Parameters"]
    process_settings.ValidateAndAssignDefaults(default_settings)
    computing_model_part = Model[process_settings["model_part_name"].GetString()]

    Logger.PrintInfo("SetCylindricalLocalAxesProcess:: ","Setting the oriented local axes...")
    process_settings.RemoveValue("model_part_name")
    return SMA.SetCylindricalLocalAxesProcess(computing_model_part, process_settings)




