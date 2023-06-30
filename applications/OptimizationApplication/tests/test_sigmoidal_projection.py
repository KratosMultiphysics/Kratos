import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.KratosUnittest import TestCase
class TestSigmoidalProjection(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()

        cls.test_model_part =  cls.model.CreateModelPart("test")
        cls.test_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 3)

        cls.test_model_part.CreateNewNode(1, 0.0,0.0,0.0)

    def _GetContainerExpression(self):
        return Kratos.Expression.NodalExpression(self.test_model_part)

    def test_ProjectForward(self):
        test_field = self._GetContainerExpression()
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1],[2],25.0,1)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1,2],[3,3],25.0,1)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1,1],[2,3],25.0,1)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1,2],[1],25.0,1)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1],[3,4],25.0,1)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[2,1],[3,4],25.0,1)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1,2],[4,3],25.0,1)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 1.5)
        forward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1,2],[3,4],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(forward_projected_field), 3.5, 4)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 5.0)
        forward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1,2],[3,4],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(forward_projected_field), 4.0, 4)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, -10.0)
        forward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1,2],[3,4],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(forward_projected_field), 3.0, 4)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 1.51)
        forward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1,2],[3,4],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(forward_projected_field), 3.6224593312018545, 10)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 1.49)
        forward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1,2],[3,4],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(forward_projected_field), 3.3775406687981455, 10)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 2.5)
        forward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1,2,3],[4,5,6],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(forward_projected_field), 5.5, 4)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 2.0)
        forward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1,2,3],[4,5,6],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(forward_projected_field), 5.0, 4)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 1.0)
        forward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1,2,3],[4,5,6],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(forward_projected_field), 4.0, 4)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 3.0)
        forward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(test_field,[1,2,3],[4,5,6],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(forward_projected_field), 6.0, 4)

    def test_ProjectBackward(self):
        test_field = self._GetContainerExpression()
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1],[2],25.0,1)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1,2],[3,3],25.0,1)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1,1],[2,3],25.0,1)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1,2],[1],25.0,1)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1],[3,4],25.0,1)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[2,1],[3,4],25.0,1)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1,2],[4,3],25.0,1)
        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 5.5)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1,2],[3,4],25.0,1)
        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 2.5)
        with self.assertRaises(RuntimeError):
            _ = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1,2],[3,4],25.0,1)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 5.5)
        backward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1,2,3],[4,5,6],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(backward_projected_field), 2.5, 4)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 5.999999999)
        backward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1,2,3],[4,5,6],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(backward_projected_field), 2.914465315084121, 10)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 4.5)
        backward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1,2,3],[4,5,6],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(backward_projected_field), 1.5, 4)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 4.0)
        backward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1,2,3],[4,5,6],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(backward_projected_field), 1.0, 4)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 5.0)
        backward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1,2,3],[4,5,6],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(backward_projected_field), 2.0, 4)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 800.0)
        backward_projected_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(test_field,[1,2,3],[400,420,1000],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(backward_projected_field), 2.512837077723448, 10)

    def test_ComputeFirstDerivative(self):
        test_field = self._GetContainerExpression()
        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 5.5)
        derivative_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.CalculateForwardProjectionGradient(test_field,[1,2,3],[4,5,6],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(derivative_field), 0.0, 4)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 2.5)
        derivative_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.CalculateForwardProjectionGradient(test_field,[1,2,3],[4,5,6],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(derivative_field), 12.5, 4)

        Kratos.Expression.LiteralExpressionIO.SetData(test_field, 1.5)
        derivative_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.CalculateForwardProjectionGradient(test_field,[1,2,3],[0,500,600],25.0,1)
        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(derivative_field), 6250.0, 4)

if __name__ == '__main__':
    Kratos.KratosUnittest.main()