import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam

## In this case, the scalar value is automatically fixed.

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeNodalYoungModulusProcess(Model, settings["Parameters"])

class ImposeNodalYoungModulusProcess(Process):

    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]
        settings.AddEmptyValue("constrained").SetBool(True)

        self.process = KratosDam.DamNodalYoungModulusProcess(model_part, settings)


    def ExecuteBeforeSolutionLoop(self):

        self.process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
