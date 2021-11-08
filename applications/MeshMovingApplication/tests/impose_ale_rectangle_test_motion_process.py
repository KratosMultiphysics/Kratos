import KratosMultiphysics as KM
from math import cos, sin

def Factory(settings, Model):
    if(not isinstance(settings, KM.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeALERectangleTestMotionProcess(Model, settings["Parameters"])

class ImposeALERectangleTestMotionProcess(KM.Process):
    def __init__(self, Model, settings):
        KM.Process.__init__(self)
        self.model_part = Model[settings["model_part_name"].GetString()]

    def ExecuteInitialize(self):
        KM.VariableUtils().ApplyFixity(KM.MESH_DISPLACEMENT_X, True, self.model_part.Nodes)
        KM.VariableUtils().ApplyFixity(KM.MESH_DISPLACEMENT_Y, True, self.model_part.Nodes)
        KM.VariableUtils().ApplyFixity(KM.MESH_DISPLACEMENT_Z, True, self.model_part.Nodes)

    def ExecuteInitializeSolutionStep(self):
        xc = 0.375
        yc = 0.250
        time = self.model_part.ProcessInfo[KM.TIME]
        pi = 3.141592653589793
        ut = 0.10 * sin(pi * time)
        phi = 0.10 * sin(2.0 * pi * time)
        c = cos(phi)
        s = sin(phi)
        for node in self.model_part.Nodes:
            rx = node.X0 - xc
            ry = node.Y0 - yc
            ur = xc + c * rx - s * ry - node.X0
            vr = yc + s * rx + c * ry - node.Y0
            node.SetSolutionStepValue(KM.MESH_DISPLACEMENT_X, ut + ur)
            node.SetSolutionStepValue(KM.MESH_DISPLACEMENT_Y, vr)

