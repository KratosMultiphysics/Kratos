import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.KratosUnittest as UnitTest

class TestSmoothClamper(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.min = -150
        cls.max = 150

        for i in range(101):
            node: Kratos.Node = cls.model_part.CreateNewNode(i + 1, i + 1, i + 2, i + 3)
            v = cls.min + float(i) * (cls.max - cls.min)  / 100.0
            node.SetValue(Kratos.PRESSURE, v)

        cls.x_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(cls.model_part.Nodes, Kratos.PRESSURE)
        cls.x_ta.Check()
        cls.x_ta.CollectData()
        cls.clamper = KratosSI.SmoothClamper(-10, 10)

    def test_Clamp(self) -> None:
        y = self.clamper.ProjectForward(self.x_ta)
        for i, y_i in enumerate(y.data):
            if self.x_ta.data[i] < -10.0:
                self.assertEqual(y_i, -10.0)
            elif self.x_ta.data[i] > 10.0:
                self.assertEqual(y_i, 10.0)
        self.assertAlmostEqual(np.linalg.norm(y.data), 98.78158380993898)

    def test_InverseClamp(self) -> None:
        y = self.clamper.ProjectForward(self.x_ta)
        x = self.clamper.ProjectBackward(y)

        x_clamped = np.clip(self.x_ta.data, -10, 10)
        self.assertAlmostEqual(np.linalg.norm(x.data - x_clamped), 0.0)

    def test_ClampDerivative(self) -> None:
        delta = 1e-8
        y_ref = self.clamper.ProjectForward(self.x_ta)

        perturbed_x_ta = self.x_ta.Clone()
        perturbed_x_ta.data[:] += delta
        y = self.clamper.ProjectForward(perturbed_x_ta)
        dy_dx = self.clamper.CalculateForwardProjectionGradient(self.x_ta)

        for i, dy_dx in enumerate(dy_dx.data):
            self.assertAlmostEqual(dy_dx, (y.data[i] - y_ref.data[i]) / delta, 5)

if __name__ == '__main__':
    UnitTest.main()