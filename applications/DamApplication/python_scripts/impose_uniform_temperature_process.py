from KratosMultiphysics import *
from KratosMultiphysics.PoromechanicsApplication import *

## This proces sets the value of a scalar variable using the ApplyConstantScalarValueProcess.
## In this case, the scalar value is not automatically fixed, so the fixicity must be introduced before.

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeUniformTemperatureProcess(Model, settings["Parameters"])

class ImposeUniformTemperatureProcess(Process):
    
    def __init__(self, Model, settings ):

        Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]
        
        params = Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("mesh_id",settings["mesh_id"])
        params.AddValue("is_fixed",settings["is_fixed"])
        params.AddValue("variable_name",settings["variable_name"])
        params.AddValue("value",settings["value"])

        if settings["table"].GetInt() == 0:
            self.process = ApplyConstantScalarValueProcess(model_part, params)
        else:
            params.AddValue("table",settings["table"])
            self.process = ApplyDoubleTableProcess(model_part, params)
                 

    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
