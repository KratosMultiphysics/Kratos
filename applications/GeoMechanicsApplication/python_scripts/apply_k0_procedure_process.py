import KratosMultiphysics as Core
import KratosMultiphysics.GeoMechanicsApplication as Geo


def Factory(settings, model):
    if not isinstance(settings, Core.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")

    model_part_name = settings["Parameters"]["model_part_name"]
    model_part = model[model_part_name.GetString()]
    params = Core.Parameters("{}")
    params.AddValue("model_part_name", model_part_name)
    return Geo.ApplyK0ProcedureProcess(model_part, params)
