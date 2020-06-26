from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *


def Factory(settings, Model):
    if not isinstance(settings, Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeTemperaturebyDeviceProcess(Model, settings["Parameters"])

class ImposeTemperaturebyDeviceProcess(Process):

    def __init__(self, Model, settings ):

        Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]

        self.process = DamTemperaturebyDeviceProcess(model_part, settings)

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
