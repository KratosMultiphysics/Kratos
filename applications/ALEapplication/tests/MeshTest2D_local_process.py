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

        self.movement_type = Parameters["movement_type"].GetString()
        self.amplification_factor = Parameters["amplification_factor"].GetDouble()
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
              node.Fix(KratosALE.MESH_DISPLACEMENT_Z)
            if(node.Y > 0.4999 or node.Y < 0.0001):
              node.Fix(KratosALE.MESH_DISPLACEMENT_Y)
              node.Fix(KratosALE.MESH_DISPLACEMENT_Z)


    def ExecuteBeforeSolutionLoop(self):
        pass


    def ExecuteInitializeSolutionStep(self):
        time = self.interface_model_part.ProcessInfo[KratosMultiphysics.TIME]
        factor = time*self.amplification_factor

        for node in self.interface_model_part.Nodes:
            # Compute the prescribed values
            if self.movement_type == "translation_x":
                valueX = 0.5 * factor
                valueY = 0
            elif self.movement_type == "translation_y":
                valueX = 0
                valueY = 0.2 * factor
            elif self.movement_type == "bending":
                valueX = 0
                valueY = 3.0 * node.X0 * node.X0 * factor
            elif self.movement_type == "rotation":
                xOld = node.X0
                yOld = node.Y0
                valueX = (math.cos(factor*.6)*xOld + math.sin(factor*.6)*yOld) - xOld
                valueY = (-math.sin(factor*.6)*xOld + math.cos(factor*.6)*yOld) - yOld
            else:
                wait = input("Wrong type of movement specified, please correct input")

            # Set the prescribed values
            node.SetSolutionStepValue(KratosALE.MESH_DISPLACEMENT_X,0,valueX)
            node.SetSolutionStepValue(KratosALE.MESH_DISPLACEMENT_Y,0,valueY)


    def ExecuteFinalizeSolutionStep(self):
        pass


    def ExecuteBeforeOutputStep(self):
        pass


    def ExecuteAfterOutputStep(self):
        pass


    def ExecuteFinalize(self):
        pass
