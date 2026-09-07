import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam

## In this case, the nodal value is assigned as a prescribed material field.
## It is not a degree of freedom, so it is not fixed.

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeNodalYoungModulusProcess(Model, settings["Parameters"])

class ImposeNodalYoungModulusProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]

        self.process = KratosDam.DamNodalYoungModulusProcess(model_part, settings)


    def ExecuteBeforeSolutionLoop(self):

        self.process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
