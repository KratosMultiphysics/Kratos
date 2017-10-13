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
    return EmbeddedSlipBoundaryConditionTestLocalProcess(Model, settings["Parameters"])

class EmbeddedSlipBoundaryConditionTestLocalProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):

    def __init__(self,model_part,params):

        self.parameters = params
        self.fluid_model_part = model_part["MainModelPart"]
        self.distance_inlet = params["distance_inlet"].GetDouble()
        self.distance_outlet = params["distance_outlet"].GetDouble()


    def ExecuteInitialize(self):

        # Set the distance function
        for node in self.fluid_model_part.Nodes:
            distance = node.Z - ((self.distance_outlet-self.distance_inlet)*node.X*node.X/25.0 + self.distance_inlet)
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, distance)

        # Deactivate the elements that have negative distance value
        for elem in self.fluid_model_part.Elements:
            interior = 0
            for node in elem.GetNodes():
                if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0:
                    interior += 1
            if interior == 4:
                elem.Set(KratosMultiphysics.ACTIVE,False)

        # Set the top velocity constraint
        for node in self.fluid_model_part.GetSubModelPart("Top3D").Nodes:
            vel = KratosMultiphysics.Vector(3)
            vel[0] = 0.0
            vel[1] = 0.0
            vel[2] = 0.0

            node.Fix(KratosMultiphysics.VELOCITY_Z)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, vel)

        # Set the inlet
        for node in self.fluid_model_part.GetSubModelPart("Inlet3D").Nodes:
            if (node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0.0):
                vel = KratosMultiphysics.Vector(3)
                vel[0] = 1.0
                vel[1] = 0.0
                vel[2] = 0.0

                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.Fix(KratosMultiphysics.VELOCITY_Z)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, vel)

        # Set the SLIP flag
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, True, self.fluid_model_part.Elements)


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
        pass
