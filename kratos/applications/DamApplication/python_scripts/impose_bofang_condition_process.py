from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

## This proces sets the value of a scalar variable using the BofangConditionTemperatureProcess.
## In this case, the scalar value is automatically fixed.

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeBofangConditionProcess(Model, settings["Parameters"])

class ImposeBofangConditionProcess(Process):
    
    def __init__(self, Model, settings ):

        Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]
        
        settings.AddEmptyValue("is_fixed").SetBool(True)
        
        # TODO: Maybe we can add a new parameter to modify in time.   
        
        self.process = BofangConditionTemperatureProcess(model_part, settings) 
                 

    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
