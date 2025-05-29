import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam

## In this case, the scalar value is automatically fixed.

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysicsParameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeChemoMechanicalAgingProcess(Model, settings["Parameters"])

class ImposeChemoMechanicalAgingProcess(KratosMultiphysicsProcess):

    def __init__(self, Model, settings ):

        KratosMultiphysicsProcess.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]
        self.process = KratosDam.DamChemoMechanicalAgingYoungProcess(model_part, settings)

    def ExecuteBeforeSolutionLoop(self):
        self.process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.process.ExecuteInitializeSolutionStep()
