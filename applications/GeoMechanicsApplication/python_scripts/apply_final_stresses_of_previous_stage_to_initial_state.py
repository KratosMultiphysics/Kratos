import KratosMultiphysics as Core
import KratosMultiphysics.GeoMechanicsApplication as Geo


def Factory(settings, model):
    if not isinstance(settings, Core.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")

    model_part = model[settings["Parameters"]["model_part_name"].GetString()]
    return Geo.ApplyFinalStressesOfPreviousStageToInitialState(model_part, settings["Parameters"])