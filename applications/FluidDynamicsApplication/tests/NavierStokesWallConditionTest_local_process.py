from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

KratosMultiphysics.CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLocalProcess(Model, settings["Parameters"])

class ApplyLocalProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):

    def __init__(self,model_part,params):

        self.ext_pres_res = params["ext_pres_res"].GetDouble()
        self.fluid_model_part = model_part[params["fluid_model_part_name"].GetString()]


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

        react_x = 0.0

        for node in self.fluid_model_part.Nodes:
            react_x += node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)

        self.assertAlmostEqual(-react_x, self.ext_pres_res, delta=5e-3)
