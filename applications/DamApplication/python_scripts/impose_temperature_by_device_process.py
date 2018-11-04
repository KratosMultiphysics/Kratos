from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

## This proces sets the value of a scalar variable using the BofangConditionTemperatureProcess.
## In this case, the scalar value is automatically fixed.

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
