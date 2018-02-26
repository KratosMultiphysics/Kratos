from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

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
 
        if settings["table"].GetInt() == 0:
            param = Parameters("{}")
            param.AddValue("model_part_name",settings["model_part_name"])
            param.AddValue("mesh_id",settings["mesh_id"])
            param.AddEmptyValue("is_fixed").SetBool(False)
            param.AddValue("variable_name",settings["variable_name"])
            param.AddValue("value",settings["value"])
            self.process = ApplyConstantScalarValueProcess(model_part, param) 
        else:
            self.process = DamFixTemperatureConditionProcess(model_part, settings)

    def ExecuteInitialize(self):
        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.process.ExecuteInitializeSolutionStep()
            
    def ExecuteFinalizeSolutionStep(self):
        self.process.ExecuteFinalizeSolutionStep()
