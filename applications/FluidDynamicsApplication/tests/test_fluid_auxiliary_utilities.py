import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtils

class FluidAuxiliaryUtilitiesTest(UnitTest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        fluid_model_part = self.model.CreateModelPart("FluidModelPart")
        fluid_model_part.SetBufferSize(2)
        fluid_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 2)
        fluid_model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        fluid_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        Kratos.ModelPartIO("Cavity/square5").ReadModelPart(fluid_model_part)

    def testCalculateFluidVolume(self):
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        fluid_volume = KratosFluid.FluidAuxiliaryUtilities.CalculateFluidVolume(fluid_model_part)
        self.assertAlmostEqual(fluid_volume, 1.0, 12)

    def testCalculateFluidPositiveVolume(self):
        # Set fluid level set
        level_set_y = 1.0/3.0
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        for node in fluid_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISTANCE, 0, node.Y - level_set_y)

        # Calculate the fluid positive volume
        fluid_positive_volume = KratosFluid.FluidAuxiliaryUtilities.CalculateFluidPositiveVolume(fluid_model_part)
        self.assertAlmostEqual(fluid_positive_volume, 1.0 - level_set_y, 12)

    def testCalculateFluidNegativeVolume(self):
        # Set fluid level set
        level_set_y = 2.0/3.0
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        for node in fluid_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISTANCE, 0, node.Y - level_set_y)

        # Calculate the fluid negative volume
        fluid_negative_volume = KratosFluid.FluidAuxiliaryUtilities.CalculateFluidNegativeVolume(fluid_model_part)
        self.assertAlmostEqual(fluid_negative_volume, level_set_y, 12)

    def tearDown(self):
        KratosUtils.DeleteFileIfExisting("Cavity/square5.time")


if __name__ == '__main__':
    UnitTest.main()
