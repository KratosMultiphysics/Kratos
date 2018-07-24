from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

## This proces sets the value of a scalar variable using the BofangConditionTemperatureProcess.
## In this case, the scalar value is automatically fixed.

def Factory(settings, Model):
    if not isinstance(settings, Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeReservoirTemperatureConditionProcess(Model, settings["Parameters"])

class ImposeReservoirTemperatureConditionProcess(Process):

    def __init__(self, Model, settings ):

        Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]

        settings.AddEmptyValue("is_fixed").SetBool(True)

        self.components_process_list = []

        if "BOFANGTEMPERATURE" in settings["model_part_name"].GetString():
            self.components_process_list.append(DamBofangConditionTemperatureProcess(model_part, settings))

        if "CONSTANTRESERVOIRTEMPERATURE" in settings["model_part_name"].GetString():
            self.components_process_list.append(DamReservoirConstantTemperatureProcess(model_part, settings))


    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteFinalizeSolutionStep()
