from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

def Factory(settings, Model):
    if(not isinstance(settings,Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeDamSurfaceNodeProcess(Model, settings["Parameters"])

class ImposeDamSurfaceNodeProcess(Process):
    def __init__(self, Model, settings ):

        Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]

        self.process = DamSurfaceNodeProcess(model_part, settings)

    def ExecuteInitialize(self):
        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        pass
