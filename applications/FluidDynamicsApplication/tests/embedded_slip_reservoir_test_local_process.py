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
    return EmbeddedSlipReservoirTestLocalProcess(Model, settings["Parameters"])

class EmbeddedSlipReservoirTestLocalProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):

    def __init__(self,model_part,params):

        self.fluid_model_part = model_part[params["fluid_model_part_name"].GetString()]
        self.distance = params["distance"].GetDouble()


    def ExecuteInitialize(self):

        # Set the distance function
        for node in self.fluid_model_part.Nodes:
            distance = node.Z - self.distance + 8.0*node.X*node.Y*(node.X-1.0)*(node.Y-1.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, distance)

        # Deactivate the elements that have negative distance value
        for elem in self.fluid_model_part.Elements:
            interior = 0
            for node in elem.GetNodes():
                if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0:
                    interior += 1
            if interior == 4:
                elem.Set(KratosMultiphysics.ACTIVE,False)

        # Flag the elements as SLIP
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

        for node in self.fluid_model_part.Nodes:
            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0.0:

                density = node.GetSolutionStepValue(KratosMultiphysics.DENSITY,0)
                gravity = node.GetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Z,0)

                expected_solution = (2-node.Z)*density*abs(gravity)
                obtained_solution = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE,0)
                self.assertAlmostEqual(obtained_solution,expected_solution,delta=1e-01)

                obtained_solution = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0)
                self.assertAlmostEqual(obtained_solution,0.0,delta=2.5e-02)

                obtained_solution = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0)
                self.assertAlmostEqual(obtained_solution,0.0,delta=2.5e-02)

                obtained_solution = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0)
                self.assertAlmostEqual(obtained_solution,0.0,delta=2.5e-02)
