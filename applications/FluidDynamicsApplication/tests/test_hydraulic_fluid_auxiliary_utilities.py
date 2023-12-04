import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtils

class HydraulicFluidAuxiliaryUtilitiesTest(UnitTest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        fluid_model_part = self.model.CreateModelPart("FluidModelPart")
        fluid_model_part.SetBufferSize(2)
        fluid_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 3)
        fluid_model_part.AddNodalSolutionStepVariable(Kratos.AUX_DISTANCE)
        # Kratos.ModelPartIO("two_fluid_hydraulic_solver_test/3D_geometry").ReadModelPart(fluid_model_part)
        Kratos.ModelPartIO("two_fluid_hydraulic_solver_test/geo3d").ReadModelPart(fluid_model_part)


    def testCalculateWettedPerimeter(self):
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        inlet_skin_model_part = self.model.GetModelPart("FluidModelPart.AutomaticInlet3D_Inlet")
        level_set_y= 0.3

        process = Kratos.FindGlobalNodalNeighboursProcess(fluid_model_part.GetCommunicator().GetDataCommunicator(), fluid_model_part)
        process.Execute()

        node_id_map = process.GetNeighbourIds(fluid_model_part.Nodes)
        for node in inlet_skin_model_part.Nodes:
            node.Set(Kratos.INLET, True)
        for condition in inlet_skin_model_part.Conditions:
            condition.Set(Kratos.INLET, True)
        for node in inlet_skin_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.AUX_DISTANCE, 0, node.Y - level_set_y)
        WettedPerimeter= KratosFluid.HydraulicFluidAuxiliaryUtilities.CalculateWettedPetimeter(fluid_model_part,Kratos.INLET)
        theoretical_wetted_perimeter= 1.0 +2*(level_set_y)
        self.assertAlmostEqual(WettedPerimeter, theoretical_wetted_perimeter, 12)

if __name__ == '__main__':
    UnitTest.main()