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
        
        settings.AddEmptyValue("is_fixed").SetBool(True)
        
        self.components_process_list = []

        self.components_process_list.append(DamFixTemperatureConditionProcess(model_part, settings))

    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()
            
    def ExecuteFinalizeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteFinalizeSolutionStep()
