import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

## This proces sets the value of a scalar variable using the AssignScalarVariableProcess.
## In this case, the scalar value is not automatically fixed, so the fixicity must be introduced before.

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeUniformTemperatureProcess(Model, settings["Parameters"])

class ImposeUniformTemperatureProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]

        if settings["table"].GetInt() == 0:
            param = KratosMultiphysics.Parameters("{}")
            param.AddValue("model_part_name",settings["model_part_name"])
            param.AddValue("constrained",settings["constrained"])
            param.AddValue("variable_name",settings["variable_name"])
            param.AddValue("value",settings["value"])
            if settings["constrained"].GetBool():
                self.process = KratosDam.DamFixTemperatureConditionProcess(model_part, settings)
            else:
                self.process = AssignScalarVariableProcess(Model, param)
        else:
            self.process = KratosDam.DamFixTemperatureConditionProcess(model_part, settings)

    def ExecuteBeforeSolutionLoop(self):
        self.process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.process.ExecuteFinalizeSolutionStep()
