from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import math

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as KratosALE

KratosMultiphysics.CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLocalProcess(Model, settings["Parameters"])

class ApplyLocalProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):

    def __init__(self,Model,Parameters):

        self.fluid_model_part = Model[Parameters["fluid_model_part_name"].GetString()]
        self.interface_model_part = Model[Parameters["interface_model_part_name"].GetString()]


    def ExecuteInitialize(self):
        # Fix the interface mesh displacement
        for node in self.interface_model_part.Nodes:
            node.Fix(KratosALE.MESH_DISPLACEMENT_X)
            node.Fix(KratosALE.MESH_DISPLACEMENT_Y)
            node.Fix(KratosALE.MESH_DISPLACEMENT_Z)

        # Set the contour boundary conditions
        for node in self.fluid_model_part.Nodes:
            if(node.X < 0.0001 or node.X > 1.999):
              node.Fix(KratosALE.MESH_DISPLACEMENT_X)
              node.Fix(KratosALE.MESH_DISPLACEMENT_Y)
              node.Fix(KratosALE.MESH_DISPLACEMENT_Z)
            if(node.Y > 0.4999 or node.Y < 0.0001):
              node.Fix(KratosALE.MESH_DISPLACEMENT_X)
              node.Fix(KratosALE.MESH_DISPLACEMENT_Y)
              node.Fix(KratosALE.MESH_DISPLACEMENT_Z)


    def ExecuteBeforeSolutionLoop(self):
        pass


    def ExecuteInitializeSolutionStep(self):
        time = self.interface_model_part.ProcessInfo[KratosMultiphysics.TIME]

        for node in self.interface_model_part.Nodes:
            disp = KratosMultiphysics.Vector(3)
            disp[0] = 0.0    #0.05*math.sin((time-0.02)*100)
            disp[1] = 0.0
            disp[2] = 0.03*math.sin((time-0.02)*100)
            node.SetSolutionStepValue(KratosALE.MESH_DISPLACEMENT, 0, disp)

    def ExecuteFinalizeSolutionStep(self):
        pass


    def ExecuteBeforeOutputStep(self):
        pass


    def ExecuteAfterOutputStep(self):
        pass


    def ExecuteFinalize(self):
        pass
