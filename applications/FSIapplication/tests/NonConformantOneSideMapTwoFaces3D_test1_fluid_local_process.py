from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
from KratosMultiphysics import *
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

        self.positive_interface_model_part = Model[params["positive_fluid_interface_modelpart_name"].GetString()]
        self.negative_interface_model_part = Model[params["negative_fluid_interfacemodelpart_name"].GetString()]


    def ExecuteInitialize(self):
        # Assign the PRESSURE values
        for node_pos, node_neg in zip(self.positive_interface_model_part.Nodes, self.negative_interface_model_part.Nodes):
            if (node_neg.Id == node_pos.Id):
                node_pos.SetSolutionStepValue(PRESSURE, 0, 0.0) # If the node is shared, set a unique value
            else:
                node_pos.SetSolutionStepValue(PRESSURE, 0, -(-36.07*node_pos.Y*node_pos.Y+36.04*node_pos.Y-8)*(8*node_pos.Z-16*node_pos.Z*node_pos.Z)) # Positive face nodal values
                node_neg.SetSolutionStepValue(PRESSURE, 0, (-36.07*node_neg.Y*node_neg.Y+36.04*node_neg.Y-8)*(8*node_neg.Z-16*node_neg.Z*node_neg.Z))  # Negative face nodal values


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
        # Mapped VELOCITY check in positive interface
        for node in self.positive_interface_model_part.Nodes:
            obtained_velocity_value = node.GetSolutionStepValue(VELOCITY,0)
            expected_velocity_value = [-(2*node.Z-node.Y*node.Z/2.0) , 2*node.Z-node.Y*node.Z/2.0, node.X*node.Y]
            for i in range(0,3):
                self.assertAlmostEqual(obtained_velocity_value[i], expected_velocity_value[i], delta=5e-3)

        # Mapped VELOCITY check in negative interface
        for node in self.negative_interface_model_part.Nodes:
            obtained_velocity_value = node.GetSolutionStepValue(VELOCITY,0)
            expected_velocity_value = [-(2*node.Z-node.Y*node.Z/2.0) , 2*node.Z-node.Y*node.Z/2.0, node.X*node.Y]
            for i in range(0,3):
                self.assertAlmostEqual(obtained_velocity_value[i], expected_velocity_value[i], delta=5e-3)
