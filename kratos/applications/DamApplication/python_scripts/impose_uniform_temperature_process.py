from KratosMultiphysics import *

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
        
        self.process = ApplyConstantScalarValueProcess(model_part, settings) 
                 

    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
