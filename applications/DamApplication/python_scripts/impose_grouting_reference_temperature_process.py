from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

def Factory(settings, Model):
    if(not isinstance(settings,Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeGroutingReferenceTemperatureProcess(Model, settings["Parameters"])

class ImposeGroutingReferenceTemperatureProcess(Process):
    def __init__(self, Model, settings ):

        Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()
        initial_value = settings["initial_value"].GetDouble()
        time_grouting = settings["time_grouting"].GetDouble()

        self.process = DamGroutingReferenceTemperatureProcess(model_part, settings)


    def ExecuteInitialize(self):
        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.process.ExecuteInitializeSolutionStep()
