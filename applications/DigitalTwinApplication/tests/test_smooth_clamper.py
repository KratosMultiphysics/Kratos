import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.KratosUnittest as UnitTest

class TestSmoothClamper(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.min = -10.0
        cls.max = 10.0

        for i in range(101):
            node = cls.model_part.CreateNewNode(i + 1, i + 1, i + 2, i + 3)
            v = cls.min + float(i) * (cls.max - cls.min)  / 100.0
            node.SetValue(Kratos.PRESSURE, v)

        cls.x_exp = Kratos.Expression.NodalExpression(cls.model_part)
        Kratos.Expression.VariableExpressionIO.Read(cls.x_exp, Kratos.PRESSURE, False)
        cls.x = cls.x_exp.Evaluate()
        cls.clamper = KratosDT.NodeSmoothClamper(-3, 3)

    def test_Clamp(self) -> None:
        y = self.clamper.Clamp(self.x_exp).Evaluate()
        for i, y_i in enumerate(y):
            if self.x[i] < self.min:
                self.assertEqual(y_i, 0.0)
            elif self.x[i] > self.max:
                self.assertEqual(y_i, 1.0)
        self.assertAlmostEqual(np.linalg.norm(y), 6.829557700271981)

    def test_InverseClamp(self) -> None:
        y = self.clamper.Clamp(self.x_exp)
        x = self.clamper.InverseClamp(y).Evaluate()

        x_clamped = np.clip(self.x_exp.Evaluate(), -3, 3)
        self.assertAlmostEqual(np.linalg.norm(x - x_clamped), 0.0)

    def test_ClampDerivative(self) -> None:
        delta = 1e-8
        y_ref = self.clamper.Clamp(self.x_exp).Evaluate()
        y = self.clamper.Clamp(self.x_exp + delta).Evaluate()
        dy_dx = self.clamper.ClampDerivative(self.x_exp).Evaluate()

        for i, dy_dx in enumerate(dy_dx):
            self.assertAlmostEqual(dy_dx, (y[i] - y_ref[i]) / delta)

if __name__ == '__main__':
    UnitTest.main()