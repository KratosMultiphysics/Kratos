import KratosMultiphysics as KM
import KratosMultiphysics.SPHApplication as SPH
from KratosMultiphysics import Logger

def Factory(settings, Model):
    
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object")
    
    default_settings = KM.Parameters(
        """{
            "model_part_name": "set_model_part_name"
        }""");

    process_settings = settings["Parameters"]
    process_settings.ValidateAndAssignDefaults(default_settings)
    computing_model_part = Model[process_settings["model_part_name"].GetString()]

    return SPH.ComputeVolumeProcess(computing_model_part, process_settings)
