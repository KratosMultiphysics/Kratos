import KratosMultiphysics
import KratosMultiphysics.FluidTransportApplication as KratosFluidTransport

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyScalarConstraintTableProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ApplyScalarConstraintTableProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("variable_name",settings["variable_name"])
        if settings.Has("is_fixed"):
            params.AddValue("is_fixed",settings["is_fixed"])
        params.AddValue("value",settings["value"])
        if settings["table"].GetInt() == 0:
            self.process = KratosMultiphysics.ApplyConstantScalarValueProcess(model_part, params)
        else:
            params.AddValue("table",settings["table"])
            self.process = KratosFluidTransport.ApplyDoubleTableProcess(model_part, params)

    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()