import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")
    model_part = model[settings["Parameters"]["model_part_name"].GetString()]
    return KratosGeo.FindNeighbourElementsOfConditionsProcess(model_part)

