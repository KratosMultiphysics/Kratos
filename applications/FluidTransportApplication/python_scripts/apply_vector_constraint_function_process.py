import KratosMultiphysics
import KratosMultiphysics.FluidTransportApplication as KratosFluidTransport
import math

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyVectorConstraintFunctionProcess(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class ApplyVectorConstraintFunctionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        self.components_process_list = []

        if settings["active"][0].GetBool() == True:
            for node in model_part.Nodes:
                # velocity = 10000*node.Y*(1-node.X*node.X)
                velocity =  2.0 * (node.Y - 0.5)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,velocity)

        if settings["active"][1].GetBool() == True:
            for node in model_part.Nodes:
                # velocity = -10000*node.X*(1-node.Y*node.Y)
                velocity =  - 2.0 * (node.X - 0.5)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,velocity)

        if settings["active"][2].GetBool() == True:
            for node in model_part.Nodes:
                velocity = 0.0
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,velocity)

    # def ExecuteInitialize(self):

    #     for component in self.components_process_list:
    #         component.ExecuteInitialize()

    # def ExecuteInitializeSolutionStep(self):

    #     for component in self.components_process_list:
    #         component.ExecuteInitializeSolutionStep()