from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.FSIApplication import *

CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLocalProcess(Model, settings["Parameters"])

class ApplyLocalProcess(Process, KratosUnittest.TestCase):

    def __init__(self,Model,params):

        self.structure_interface_model_part = Model[params["structure_interface_modelpart_name"].GetString()]


    def ExecuteInitialize(self):
        # Parabolic velocity distribution
        for node in self.structure_interface_model_part.Nodes:
            node.SetSolutionStepValue(VELOCITY_X, 0, -(2*node.Z - node.Y*node.Z/2.0))
            node.SetSolutionStepValue(VELOCITY_Y, 0, 2*node.Z - node.Y*node.Z/2.0)
            node.SetSolutionStepValue(VELOCITY_Z, 0, node.X*node.Y)


    def ExecuteBeforeSolutionLoop(self):
        pass


    def ExecuteInitializeSolutionStep(self):
        pass


    def ExecuteFinalizeSolutionStep(self):
        pass


    def ExecuteBeforeOutputStep(self):
        pass


    def ExecuteAfterOutputStep(self):
        pass


    def ExecuteFinalize(self):
        # Mapped PRESSURE check
        for node in self.structure_interface_model_part.Nodes:
            obtained_pressure_value = node.GetSolutionStepValue(POSITIVE_FACE_PRESSURE,0)
            expected_pressure_value = -(-36.07*node.Y*node.Y+36.04*node.Y-8)*(8*node.Z-16*node.Z*node.Z)
            self.assertAlmostEqual(obtained_pressure_value, expected_pressure_value, delta=2.5e-2)

        for node in self.structure_interface_model_part.Nodes:
            obtained_pressure_value = node.GetSolutionStepValue(NEGATIVE_FACE_PRESSURE,0)
            expected_pressure_value = (-36.07*node.Y*node.Y+36.04*node.Y-8)*(8*node.Z-16*node.Z*node.Z)
            self.assertAlmostEqual(obtained_pressure_value, expected_pressure_value, delta=2.5e-2)
