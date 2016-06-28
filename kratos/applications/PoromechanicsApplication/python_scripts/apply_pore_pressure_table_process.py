from KratosMultiphysics import *
from KratosMultiphysics.PoromechanicsApplication import *


def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyPorePressureTableProcess(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class ApplyPorePressureTableProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        params = Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("mesh_id",settings["mesh_id"])
        params.AddValue("variable_name",settings["variable_name"])
        params.AddValue("is_fixed",settings["is_fixed"])
        
        if settings["hydrostatic"].GetBool() == False:
            params.AddValue("value",settings["value"])
            if settings["table"].GetInt() == 0:
                self.process = ApplyConstantScalarValueProcess(model_part, params)
            else:
                params.AddValue("table",settings["table"])
                self.process = ApplyDoubleTableProcess(model_part, params)
        else:
            params.AddValue("gravity_direction",settings["gravity_direction"])
            params.AddValue("reference_coordinate",settings["reference_coordinate"])
            params.AddValue("specific_weight",settings["specific_weight"])
            if settings["table"].GetInt() == 0:
                self.process = ApplyConstantHydrostaticPressureProcess(model_part, params)
            else:
                params.AddValue("table",settings["table"])
                self.process = ApplyHydrostaticPressureTableProcess(model_part, params)

    def ExecuteInitialize(self):
        
        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        
        self.process.ExecuteInitializeSolutionStep()