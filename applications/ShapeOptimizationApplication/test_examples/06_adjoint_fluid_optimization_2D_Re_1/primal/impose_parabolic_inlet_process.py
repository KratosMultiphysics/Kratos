import KratosMultiphysics
from math import cos

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeParabolicInletProcess(Model, settings["Parameters"])

class ImposeParabolicInletProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        
        self.model_part = Model[settings["model_part_name"].GetString()]

    def ExecuteInitialize(self):
        for node in self.model_part.Nodes:
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.Fix(KratosMultiphysics.VELOCITY_Y)
            node.Fix(KratosMultiphysics.VELOCITY_Z)
        
    def ExecuteInitializeSolutionStep(self):
        Um = 1e-4
        H = 0.41
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        if time < 0.5:
            pi = 3.141592653589793
            Um = Um * 0.5 * (1.0 - cos(pi * time / 0.5))
        for node in self.model_part.Nodes:
            y = node.Y
            u = 4.0 * Um * y * (H - y) / H**2
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,u)
