
import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")

    default_settings = KratosMultiphysics.Parameters("""
            {
                "help"                    : "This process applies absorbing boundary conditions to a modelpart.",
                "model_part_name"         : "please_specify_model_part_name",
                "absorbing_factors"       : [1,1],
                "virtual_thickness"       : 1e30,
                "skip_internal_forces"    : false
            }
            """
                                                     )
    boundary_settings = settings["Parameters"]
    boundary_settings.ValidateAndAssignDefaults(default_settings)


    return SetAbsorbingBoundaryParametersProcess(model, boundary_settings)


class SetAbsorbingBoundaryParametersProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        model_part = model[settings["model_part_name"].GetString()]

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("absorbing_factors", settings["absorbing_factors"])
        params.AddValue("virtual_thickness", settings["virtual_thickness"])
        params.AddValue("skip_internal_forces", settings["skip_internal_forces"])

        self.process = KratosGeo.SetAbsorbingBoundaryParametersProcess(model_part, params)

    def ExecuteInitialize(self):
        self.process.ExecuteInitialize()
