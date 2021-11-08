import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestImporting(KratosUnittest.TestCase):

    def test_importing(self):
        """With this test we chack that the keys are not reordered after
        importing applications. Also we check that the nodal values can
        still be accessed with the variables. This would break if the
        keys are being reordered
        """
        current_model = KM.Model()

        model_part= current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        key_VELOCITY_before_import = KM.VELOCITY.Key() # test a var from the Kernel

        model_part.CreateNewNode(1,0.0,0.0,0.0)
        model_part.SetBufferSize(3)
        model_part.CloneTimeStep(1.0)
        model_part.Nodes[1].SetSolutionStepValue(KM.VELOCITY_X,0,1.0) #here i set VELOCITY_X to 1.0

        # import aux_external_import
        self._aux_func(model_part) #here i set VELOCITY_Y to 2.0
        self.assertEqual(key_VELOCITY_before_import, KM.VELOCITY.Key())

        self.assertTrue(model_part.Nodes[1].GetSolutionStepValue(KM.VELOCITY_X) == 1.0)
        self.assertTrue(model_part.Nodes[1].GetSolutionStepValue(KM.VELOCITY_Y) == 2.0)

    def _aux_func(self,model_part):
        try:
            import KratosMultiphysics.FluidDynamicsApplication
        except:
            self.skipTest("KratosMultiphysics.FluidDynamicsApplication is not available")
        model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        model_part.Nodes[1].SetSolutionStepValue(KM.VELOCITY_Y,0,2.0)

    def test_has_application(self):
        current_model = KM.Model()

        self.assertTrue(KM.KratosGlobals.Kernel.IsImported("KratosMultiphysics"))

        try:
            import KratosMultiphysics.ExternalSolversApplication
        except:
            self.skipTest("KratosMultiphysics.ExternalSolversApplication is not available")

        self.assertTrue(KM.KratosGlobals.Kernel.IsImported("ExternalSolversApplication"))

    def test_variable_keys_reordering(self):
        key_vel_before_import = KM.VELOCITY.Key() # test a var from the Kernel

        try:
            import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
        except:
            self.skipTest("KratosMultiphysics.FluidDynamicsApplication is not available")

        key_ssp_before_import = KratosCFD.SUBSCALE_PRESSURE.Key() # test an Application-Variable

        try:
            import KratosMultiphysics.StructuralMechanicsApplication
        except:
            self.skipTest("KratosMultiphysics.StructuralMechanicsApplication is not available")

        self.assertEqual(key_vel_before_import, KM.VELOCITY.Key())
        self.assertEqual(key_ssp_before_import, KratosCFD.SUBSCALE_PRESSURE.Key())



if __name__ == '__main__':
    KratosUnittest.main()
