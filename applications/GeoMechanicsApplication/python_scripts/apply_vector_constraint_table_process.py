import KratosMultiphysics as Core
import KratosMultiphysics.GeoMechanicsApplication as Geo


def Factory(settings, Model):
    if not isinstance(settings, Core.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")

    params = settings["Parameters"]
    model_part = Model[params["model_part_name"].GetString()]
    return Geo.ApplyVectorConstraintTableProcess(model_part, params)
