
import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")
    return SetAbsorbingBoundaryParametersProcess(model, settings["Parameters"])


class SetAbsorbingBoundaryParametersProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        model_part = model[settings["model_part_name"].GetString()]

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("absorbing_factors", settings["absorbing_factors"])
        params.AddValue("virtual_thickness", settings["virtual_thickness"])

        self.process = KratosGeo.SetAbsorbingBoundaryParametersProcess(model_part, params)

    def ExecuteInitialize(self):
        self.process.ExecuteInitialize()
