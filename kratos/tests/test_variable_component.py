import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestVariableComponent(KratosUnittest.TestCase):
    def test_GetSourceVariable(self):
        component_variable = KratosMultiphysics.VELOCITY_X
        source_variable = KratosMultiphysics.VELOCITY
        false_source_variable = KratosMultiphysics.DISPLACEMENT
        self.assertEqual(component_variable.GetSourceVariable(), source_variable)
        self.assertNotEqual(component_variable.GetSourceVariable(), false_source_variable)

    def test_AdjointVariable(self):
        KM = KratosMultiphysics
        self.assertTrue(KM.PRESSURE.HasAdjointVariable())
        self.assertEqual(KM.PRESSURE.GetAdjointVariable(), KM.ADJOINT_PRESSURE)
        self.assertEqual(KM.DISPLACEMENT.GetAdjointVariable(), KM.ADJOINT_DISPLACEMENT)
        self.assertEqual(KM.DISPLACEMENT_X.GetAdjointVariable(), KM.ADJOINT_DISPLACEMENT_X)
        self.assertEqual(KM.VELOCITY_Z.GetAdjointVariable(), KM.ADJOINT_VELOCITY_Z)
        self.assertEqual(KM.ROTATION.GetAdjointVariable(), KM.ADJOINT_ROTATION)

        # no adjoint to primal link and no adjoint for ACCELERATION
        self.assertFalse(KM.ADJOINT_DISPLACEMENT.HasAdjointVariable())
        self.assertFalse(KM.ACCELERATION.HasAdjointVariable())
        self.assertFalse(KM.DENSITY.HasAdjointVariable())
        with self.assertRaisesRegex(RuntimeError, "Adjoint variable for Variable \"ACCELERATION\" was not assigned"):
            KM.ACCELERATION.GetAdjointVariable()

if __name__ == '__main__':
    KratosUnittest.main()
