from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import math
import numpy as np


def Factory(Model):
    #if(type(settings) != KratosMultiphysics.Parameters):
    #    raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return RotateModePartProcess(Model)

class RotateModePartProcess(KratosMultiphysics.Process):
    """A process for to rotate full model part."""

    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]
        for node in model_part.Nodes:
            node.Fix(DISPLACEMENT)


        # detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString() ):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("the second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        self.interval = [settings["interval"][0].GetDouble(), settings["interval"][1].GetDouble()]

        if(not settings["value"][0].IsNull()):
            self.rot_x = True
            self.factor_x = 1
        else:
            self.rot_x= False
    def ExecuteInitializeSolutionStep(self, alpha):
        time = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)

        if(time > self.interval[0] and time < self.interval[1]):
            rotation = np.array([[math.cos(alpha), -math.sin(alpha), 0],[math.sin(alpha), math.cos(alpha), 0],[0, 0, 1]])
            for node in self.model_part.Nodes:
                coords = np.array([node.X0, node.Y0, node.Z0])
                new_coords = rotation.dot(coords)
                disps = new_coords-coords
                disp = Vector(3)
                disp[0] = disps[0]
                disp[1] = disps[1]
                disp[2] = disps[2]
                node.SetSolutionStepValue(DISPLACEMENT, disp)
