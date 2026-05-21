import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")

    params = settings["Parameters"]
    model_part = Model[params["model_part_name"].GetString()]
    return KratosGeo.ApplyExcavationProcess(model_part, params)

