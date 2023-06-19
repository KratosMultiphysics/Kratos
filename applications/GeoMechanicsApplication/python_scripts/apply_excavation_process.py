
import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyExcavationProcess(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class ApplyExcavationProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("variable_name",settings["variable_name"])

        if settings.Has("deactivate_soil_part"):
            params.AddValue("deactivate_soil_part",settings["deactivate_soil_part"])
            self.process = KratosGeo.ApplyExcavationProcess(model_part, params)

    def ExecuteInitialize(self):
        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.process.ExecuteInitializeSolutionStep()
