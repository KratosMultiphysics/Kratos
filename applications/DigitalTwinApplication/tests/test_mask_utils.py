import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.DigitalTwinApplication as KratosDT

class TestMaskUtils(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        n = 100
        for i in range(n):
            node: Kratos.Node = cls.model_part.CreateNewNode(i+1, i, i+1, i+2)
            node.SetValue(Kratos.PRESSURE, node.Id % 3)
            node.SetValue(Kratos.DENSITY, node.Id % 5)
            node.SetValue(Kratos.HEAT_FLUX, 2 * (node.Id % 2))

    def test_Union(self):
        mask_1 = Kratos.Expression.NodalExpression(self.model_part)
        mask_2 = Kratos.Expression.NodalExpression(self.model_part)

        Kratos.Expression.VariableExpressionIO.Read(mask_1, Kratos.PRESSURE, False)
        Kratos.Expression.VariableExpressionIO.Read(mask_2, Kratos.DENSITY, False)

        union_mask: Kratos.Expression.NodalExpression = KratosDT.MaskUtils.Union(mask_1, mask_2, 2)
        Kratos.Expression.VariableExpressionIO.Write(union_mask, Kratos.TEMPERATURE, False)
        for node in self.model_part.Nodes:
            if (node.Id % 3 >= 2 or node.Id % 5 >= 2):
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 2)
            else:
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 0)

    def test_Intersect(self):
        mask_1 = Kratos.Expression.NodalExpression(self.model_part)
        mask_2 = Kratos.Expression.NodalExpression(self.model_part)

        Kratos.Expression.VariableExpressionIO.Read(mask_1, Kratos.PRESSURE, False)
        Kratos.Expression.VariableExpressionIO.Read(mask_2, Kratos.DENSITY, False)

        union_mask: Kratos.Expression.NodalExpression = KratosDT.MaskUtils.Intersect(mask_1, mask_2, 2)
        Kratos.Expression.VariableExpressionIO.Write(union_mask, Kratos.TEMPERATURE, False)
        for node in self.model_part.Nodes:
            if (node.Id % 3 >= 2 and node.Id % 5 >= 2):
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 2)
            else:
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 0)

    def test_Substract(self):
        mask_1 = Kratos.Expression.NodalExpression(self.model_part)
        mask_2 = Kratos.Expression.NodalExpression(self.model_part)

        Kratos.Expression.VariableExpressionIO.Read(mask_1, Kratos.PRESSURE, False)
        Kratos.Expression.VariableExpressionIO.Read(mask_2, Kratos.DENSITY, False)

        union_mask: Kratos.Expression.NodalExpression = KratosDT.MaskUtils.Substract(mask_1, mask_2, 2)
        Kratos.Expression.VariableExpressionIO.Write(union_mask, Kratos.TEMPERATURE, False)
        for node in self.model_part.Nodes:
            if (node.Id % 3 >= 2 and node.Id % 5 < 2):
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 2)
            else:
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 0)

    def test_GetMaskSize(self):
       mask = Kratos.Expression.NodalExpression(self.model_part)
       Kratos.Expression.VariableExpressionIO.Read(mask, Kratos.DENSITY, False)

       self.assertEqual(KratosDT.MaskUtils.GetMaskSize(mask), 80)
       self.assertEqual(KratosDT.MaskUtils.GetMaskSize(mask, 2), 60)
       self.assertEqual(KratosDT.MaskUtils.GetMaskSize(mask, 3), 40)

    def test_GetMaskNoThreshold(self):
        values = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(values, Kratos.HEAT_FLUX, False)

        mask = KratosDT.MaskUtils.GetMask(values)
        Kratos.Expression.VariableExpressionIO.Write(mask, Kratos.TEMPERATURE, False)

        for node in self.model_part.Nodes:
            self.assertEqual(node.GetValue(Kratos.TEMPERATURE), node.Id % 2)

    def test_GetMaskWithhreshold(self):
        values = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(values, Kratos.HEAT_FLUX, False)

        mask = KratosDT.MaskUtils.GetMask(values, 1.0)
        Kratos.Expression.VariableExpressionIO.Write(mask, Kratos.TEMPERATURE, False)

        for node in self.model_part.Nodes:
            self.assertEqual(node.GetValue(Kratos.TEMPERATURE), node.Id % 2)

if __name__ == '__main__':
    UnitTest.main()
