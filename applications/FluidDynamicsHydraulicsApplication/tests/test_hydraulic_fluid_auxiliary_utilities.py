import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsHydraulicsApplication as KratosFluidHydraulics
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtils

class HydraulicFluidAuxiliaryUtilitiesTest(UnitTest.TestCase):

    def setUp(self):
        """Generates the model part to work on."""

        # Set up the domain model parts
        self.model = Kratos.Model()
        test_model_part = self.model.CreateModelPart("TestModelPart")
        inlet_model_part = test_model_part.CreateSubModelPart("Inlet")
        test_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 3)

        # Add VELOCITY variable to be used as DOF and DISTANCE
        test_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        test_model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)

        # Modelpart: The inlet condition
        # 7 - 8 - 9
        # | \ | \ |
        # 4 - 5 - 6
        # | \ | \ |
        # 1 - 2 - 3

        # Populate model part nodes
        test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        test_model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
        test_model_part.CreateNewNode(3, 0.0, 2.0, 0.0)
        test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
        test_model_part.CreateNewNode(5, 0.0, 1.0, 1.0)
        test_model_part.CreateNewNode(6, 0.0, 2.0, 1.0)
        test_model_part.CreateNewNode(7, 0.0, 0.0, 2.0)
        test_model_part.CreateNewNode(8, 0.0, 1.0, 2.0)
        test_model_part.CreateNewNode(9, 0.0, 2.0, 2.0)
        inlet_model_part.AddNodes([1, 2, 3, 4, 5, 6, 7, 8, 9])

        # Create inlet conditions
        prop_0 = test_model_part.Properties[0]
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1, 2, 4], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [2, 5, 4], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 3, [2, 3, 5], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 4, [3, 6, 5], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 5, [4, 5, 7], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 6, [5, 8, 7], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 7, [5, 6, 8], prop_0)
        inlet_model_part.CreateNewCondition("SurfaceCondition3D3N", 8, [6, 9, 8], prop_0)

    def testCalculateWettedPerimeter(self):

        inlet_skin_model_part = self.model.GetModelPart("TestModelPart.Inlet")
        level_set_z= 0.3
        for node in inlet_skin_model_part.Nodes:
            node.Set(Kratos.INLET, True)
        for condition in inlet_skin_model_part.Conditions:
            condition.Set(Kratos.INLET, True)
        for node in inlet_skin_model_part.Nodes:
            node.SetValue(Kratos.AUX_DISTANCE, node.Z - level_set_z)
        WettedPerimeter= KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.CalculateWettedPetimeter(inlet_skin_model_part,Kratos.INLET, Kratos.AUX_DISTANCE, False)
        theoretical_wetted_perimeter= 2.0 +2*(level_set_z)
        self.assertAlmostEqual(WettedPerimeter, theoretical_wetted_perimeter, 12)

    def testCalculateWettedArea(self):

        inlet_skin_model_part = self.model.GetModelPart("TestModelPart.Inlet")
        level_set_z= 0.3
        for node in inlet_skin_model_part.Nodes:
            node.Set(Kratos.INLET, True)
        for condition in inlet_skin_model_part.Conditions:
            condition.Set(Kratos.INLET, True)
        for node in inlet_skin_model_part.Nodes:
            node.SetValue(Kratos.AUX_DISTANCE, node.Z - level_set_z)
        WettedArea= KratosFluidHydraulics.HydraulicFluidAuxiliaryUtilities.CalculateWettedArea(inlet_skin_model_part,Kratos.INLET,Kratos.AUX_DISTANCE, False)
        theoretical_wetted_area= 2.0*(level_set_z)
        self.assertAlmostEqual(WettedArea, theoretical_wetted_area, 12)

if __name__ == '__main__':
    UnitTest.main()