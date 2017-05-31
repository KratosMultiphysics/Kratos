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

    def __init__(self,model_part,params):

        self.interface_model_part = model_part[params["model_part_name"].GetString()]
        self.params = params


    def ExecuteInitialize(self):
        # Parabolic velocity distribution
        for node in self.interface_model_part.Nodes:
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
        for node in self.interface_model_part.Nodes:
            obtained_pressure_value = node.GetSolutionStepValue(PRESSURE,0)
            expected_pressure_value = node.X*node.Y*node.Z
            self.assertAlmostEqual(obtained_pressure_value, expected_pressure_value, delta=5e-3)

        # Force equilibrium check (mapped POINT_LOAD)
        expected_sum_fx = 69.54492675560007
        expected_sum_fy = 187.00000000000003
        expected_sum_fz = 352.3652168853

        obtained_sum_fx = 0.0
        obtained_sum_fy = 0.0
        obtained_sum_fz = 0.0

        for node in self.interface_model_part.Nodes:
            f_solid = node.GetSolutionStepValue(POINT_LOAD)
            obtained_sum_fx += f_solid[0]
            obtained_sum_fy += f_solid[1]
            obtained_sum_fz += f_solid[2]

        self.assertAlmostEqual(obtained_sum_fx, expected_sum_fx, delta=2.5)
        self.assertAlmostEqual(obtained_sum_fy, expected_sum_fy, delta=5.0)
        self.assertAlmostEqual(obtained_sum_fz, expected_sum_fz, delta=10.0)
