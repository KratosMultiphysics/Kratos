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

        node_id_list = [100,102,112,122,146,173,219,253,297,346,397]
        vel_x_expected = [1,0.28793,0.0070363,-0.1071,-0.13454,-0.12904,-0.093136,-0.07309,-0.056608,-0.044989,-0.02823]
        vel_y_expected = [0,0.0095154,0.032294,0.039908,0.036608,0.026803,0.011931,0.0068775,0.003177,0.0020416,0.0006979]
        press_expected = [-0.025594,-0.038788,-0.034704,-0.025906,-0.015924,-0.0071712,0.0010059,0.0033377,0.0046395,0.0052515,0.0055064]
        count = 0

        for node_id in node_id_list:
            vel_x_obtained = self.fluid_model_part.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0)
            vel_y_obtained = self.fluid_model_part.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0)
            press_obtained = self.fluid_model_part.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.PRESSURE,0)

            self.assertAlmostEqual(vel_x_obtained[count], vel_x_expected[count], 6)
            self.assertAlmostEqual(vel_y_obtained[count], vel_y_expected[count], 6)
            self.assertAlmostEqual(press_obtained[count], press_expected[count], 6)
            count += 1
