import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo


def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")

    params = settings["Parameters"]
    model_part = model[params["model_part_name"].GetString()]
    return KratosGeo.GeoExtrapolateIntegrationPointValuesToNodesProcess(model_part, params)
