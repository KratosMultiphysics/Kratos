import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo


def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")

    params = settings["Parameters"]
    return KratosGeo.GeoExtrapolateIntegrationPointValuesToNodesProcess(model, params)
