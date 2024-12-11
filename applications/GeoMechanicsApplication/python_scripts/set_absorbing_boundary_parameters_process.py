
import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")
    return SetAbsorbingBoundaryParametersProcess(Model, settings["Parameters"])


class SetAbsorbingBoundaryParametersProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("absorbing_factors", settings["absorbing_factors"])
        params.AddValue("virtual_thickness", settings["virtual_thickness"])
        params.AddValue("skip_internal_forces", settings["skip_internal_forces"])

        self.process = KratosGeo.SetAbsorbingBoundaryParametersProcess(model_part, params)

    def ExecuteInitialize(self):
        self.process.ExecuteInitialize()
