from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
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
            if settings.Has("gradual_excavation"):
                if settings["gradual_excavation"]:
                    params.AddValue("gravity_direction",settings["gravity_direction"])
                    self.process = KratosGeo.ApplyGradualExcavationProcess(model_part, params)
                else:
                    self.process = KratosGeo.ApplyExcavationProcess(model_part, params)
            else:
                self.process = KratosGeo.ApplyExcavationProcess(model_part, params)

    def ExecuteInitialize(self):
        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.process.ExecuteInitializeSolutionStep()