import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam

## This process sets the value of a scalar variable using the BofangConditionTemperatureProcess.
## In this case, the scalar value is automatically fixed.

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeReservoirTemperatureConditionProcess(Model, settings["Parameters"])

class ImposeReservoirTemperatureConditionProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]

        settings.AddEmptyValue("constrained").SetBool(True)

        self.components_process_list = []

        if "BOFANGTEMPERATURE" in settings["model_part_name"].GetString():
            self.components_process_list.append(KratosDam.DamBofangConditionTemperatureProcess(model_part, settings))

        if "CONSTANTRESERVOIRTEMPERATURE" in settings["model_part_name"].GetString():
            self.components_process_list.append(KratosDam.DamReservoirConstantTemperatureProcess(model_part, settings))

        if "MONITORINGRESERVOIRTEMPERATURE" in settings["model_part_name"].GetString():
            self.components_process_list.append(KratosDam.DamReservoirMonitoringTemperatureProcess(model_part, settings))

    def ExecuteBeforeSolutionLoop(self):

        for component in self.components_process_list:
            component.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteFinalizeSolutionStep()
