import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsHydraulicsApplication as KratosFluidHydraulics

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtils
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
class HydraulicFluidAuxiliaryUtilitiesTest(UnitTest.TestCase):

    def setUp(self):
        # Read test Model Part 
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("MainModelPart")
        self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        self.model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)        
        self.model_part.AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        input_filename= "channel_test/channel3D"   
        Kratos.ModelPartIO(input_filename).ReadModelPart(self.model_part)

    def testCalculateWettedPerimeter(self):
        inlet_skin_model_part = self.model.GetModelPart("MainModelPart.AutomaticInlet3D_Inlet")
        level_set_z= 0.3
        for node in inlet_skin_model_part.Nodes:
            node.Set(Kratos.INLET, True)
        for condition in inlet_skin_model_part.Conditions:
            condition.Set(Kratos.INLET, True)
        for node in inlet_skin_model_part.Nodes:
            node.SetValue(Kratos.AUX_DISTANCE, node.Z - level_set_z)
        WettedPerimeter= KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.CalculateWettedPetimeter(inlet_skin_model_part,Kratos.INLET, Kratos.AUX_DISTANCE, False)
        theoretical_wetted_perimeter= 0.5 +2*(level_set_z)
        self.assertAlmostEqual(WettedPerimeter, theoretical_wetted_perimeter, 12)

    def testCalculateWettedArea(self):
        inlet_skin_model_part = self.model.GetModelPart("MainModelPart.AutomaticInlet3D_Inlet")
        level_set_z= 0.3
        for node in inlet_skin_model_part.Nodes:
            node.Set(Kratos.INLET, True)
        for condition in inlet_skin_model_part.Conditions:
            condition.Set(Kratos.INLET, True)
        for node in inlet_skin_model_part.Nodes:
            node.SetValue(Kratos.AUX_DISTANCE, node.Z - level_set_z)
        WettedArea= KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.CalculateWettedArea(inlet_skin_model_part,Kratos.INLET,Kratos.AUX_DISTANCE, False)
        theoretical_wetted_area= 0.5*(level_set_z)
        self.assertAlmostEqual(WettedArea, theoretical_wetted_area, 12)
  
    def testFixCornerVelocity(self):
        corner_angle=10
        Kratos.NormalCalculationUtils().CalculateOnSimplex(self.model_part,3)
        Kratos.FindConditionsNeighboursProcess(self.model_part,3,3).Execute()
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X,self.model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y,self.model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Z,self.model_part)
        KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.FixCornerNodeVelocity(self.model_part , corner_angle)
        tested_node = self.model_part.GetNode(33)
        self.assertAlmostEqual(tested_node.GetSolutionStepValue(Kratos.BODY_FORCE_Z),0.0,12)

    def testNonAirGravity(self):
        for node in self.model_part.Nodes:
            distance= node.Z-0.2
            node.SetSolutionStepValue(Kratos.DISTANCE,distance)
            node.SetSolutionStepValue(Kratos.BODY_FORCE, [0.0, 0.0, -9.81])
        KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.TurnOffGravityOnAirElements(self.model_part)
        tested_node = self.model_part.GetNode(33)
        self.assertAlmostEqual(tested_node.GetSolutionStepValue(Kratos.BODY_FORCE_Z),0.0,12)


    def testMaximumWaterDepthChange(self):
        inlet_skin_model_part = self.model.GetModelPart("MainModelPart.AutomaticInlet3D_Inlet")
        for node in inlet_skin_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISTANCE, node.Z-0.3)
            node.SetValue(Kratos.AUX_DISTANCE, node.Z-0.1)
        max_depth_change= KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.MaximumWaterDepthChange(self.model_part)
        self.assertTrue(max_depth_change)
        

if __name__ == '__main__':
    UnitTest.main()