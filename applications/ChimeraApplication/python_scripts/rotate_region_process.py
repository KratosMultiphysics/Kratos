import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
    return KratosChimera.RotateRegionProcess(model_part, settings["Parameters"])
