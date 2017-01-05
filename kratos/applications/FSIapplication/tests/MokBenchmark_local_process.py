from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.FSIApplication import *
from KratosMultiphysics.ALEApplication import *

CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLocalProcess(Model, settings["Parameters"])

class ApplyLocalProcess(Process, KratosUnittest.TestCase):

    def __init__(self,Model,settings):
        self.fluid_wall_model_part = Model[settings["model_part_name"].GetString()]

    def ExecuteInitialize(self):
        pass

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
        # Get the obtained results at the control point A (top left wall corner)
        node_A = self.fluid_wall_model_part.Nodes[2122]
        pres_A = node_A.GetSolutionStepValue(PRESSURE, 0)
        disX_A = node_A.GetSolutionStepValue(MESH_DISPLACEMENT_X, 0)
        disY_A = node_A.GetSolutionStepValue(MESH_DISPLACEMENT_Y, 0)

        # Get the obtained results at the control point B (center left wall point)
        node_B = self.fluid_wall_model_part.Nodes[3100]
        pres_B = node_B.GetSolutionStepValue(PRESSURE, 0)
        disX_B = node_B.GetSolutionStepValue(MESH_DISPLACEMENT_X, 0)
        disY_B = node_B.GetSolutionStepValue(MESH_DISPLACEMENT_Y, 0)

        # Compare the control point A results against the expected ones
        self.assertAlmostEqual(pres_A, 2.5543, delta=0.05)
        self.assertAlmostEqual(disX_A, 0.062127, delta=0.0001)
        self.assertAlmostEqual(disY_A, -0.0081185, delta=0.0001)

        # Compare the control point B results against the expected ones
        self.assertAlmostEqual(pres_B, 6.7629, delta=0.05)
        self.assertAlmostEqual(disX_B, 0.022468, delta=0.0001)
        self.assertAlmostEqual(disY_B, -0.0017516, delta=0.0001)
