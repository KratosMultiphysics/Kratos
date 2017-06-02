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
    return EmbeddedArtificialCompressibilityTestLocalProcess(Model, settings["Parameters"])

class EmbeddedArtificialCompressibilityTestLocalProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):

    def __init__(self,model_part,params):

        self.fluid_model_part = model_part[params["model_part_name"].GetString()]


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

        node_id_list = [100,102,112,122,146,173,219,253,297,346,397]
        vel_x_expected = [1,0.38577,0.038665,-0.093095,-0.11588,-0.10077,-0.069464,-0.0541,-0.043088,-0.034141,-0.022389]
        vel_y_expected = [0,0.014781,0.035381,0.038773,0.0295,0.018162,0.0069178,0.0034694,0.0017436,0.00083193,0.00027465]
        press_expected = [-0.024858,-0.03033,-0.029525,-0.022382,-0.013488,-0.0066768,-0.00093377,0.0006992,0.001516,0.0019067,0.002047]
        count = 0

        for node_id in node_id_list:
            vel_x_obtained = self.fluid_model_part.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0)
            vel_y_obtained = self.fluid_model_part.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0)
            press_obtained = self.fluid_model_part.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.PRESSURE,0)

            self.assertAlmostEqual(vel_x_obtained, vel_x_expected[count], 4)
            self.assertAlmostEqual(vel_y_obtained, vel_y_expected[count], 4)
            self.assertAlmostEqual(press_obtained, press_expected[count], 4)
            count += 1
