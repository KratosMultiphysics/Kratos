import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.KratosUnittest as UnitTest

class TestSmoothClamper(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.min = -1.5
        cls.max = 1.5

        for i in range(101):
            node = cls.model_part.CreateNewNode(i + 1, i + 1, i + 2, i + 3)
            v = cls.min + float(i) * (cls.max - cls.min)  / 100.0
            node.SetValue(Kratos.PRESSURE, v)

        cls.x_exp = Kratos.Expression.NodalExpression(cls.model_part)
        Kratos.Expression.VariableExpressionIO.Read(cls.x_exp, Kratos.PRESSURE, False)
        cls.x = cls.x_exp.Evaluate()
        cls.clamper = KratosSI.NodeSmoothClamper(-10, 10)

    def test_Clamp(self) -> None:
        y = self.clamper.ProjectForward(self.x_exp).Evaluate()
        for i, y_i in enumerate(y):
            if self.x[i] < 0.0:
                self.assertEqual(y_i, -10.0)
            elif self.x[i] > 1.0:
                self.assertEqual(y_i, 10.0)
        self.assertAlmostEqual(np.linalg.norm(y), 91.57361915207677)

    def test_InverseClamp(self) -> None:
        y = self.clamper.ProjectForward(self.x_exp)
        x = self.clamper.ProjectBackward(y).Evaluate()

        x_clamped = np.clip(self.x_exp.Evaluate(), 0, 1)
        self.assertAlmostEqual(np.linalg.norm(x - x_clamped), 0.0)

    def test_ClampDerivative(self) -> None:
        delta = 1e-8
        y_ref = self.clamper.ProjectForward(self.x_exp).Evaluate()
        y = self.clamper.ProjectForward(self.x_exp + delta).Evaluate()
        dy_dx = self.clamper.CalculateForwardProjectionGradient(self.x_exp).Evaluate()

        for i, dy_dx in enumerate(dy_dx):
            self.assertAlmostEqual(dy_dx, (y[i] - y_ref[i]) / delta, 5)

if __name__ == '__main__':
    UnitTest.main()