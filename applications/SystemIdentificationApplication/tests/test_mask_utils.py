import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import ExpressionUnionType

class TestMaskUtils(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        cls.n = 10
        for i in range(cls.n):
            node: Kratos.Node = cls.model_part.CreateNewNode(i+1, i, i+1, i+2)
            node.SetValue(Kratos.PRESSURE, node.Id % 3)
            node.SetValue(Kratos.DENSITY, node.Id % 5)
            node.SetValue(Kratos.HEAT_FLUX, 2 * (node.Id % 2))

        prop = cls.model_part.CreateNewProperties(1)
        for i in range(cls.n-1):
            cls.model_part.CreateNewElement("Element3D2N", i+1, [i+1, i+2], prop).SetValue(Kratos.DENSITY, i % 2)
            cls.model_part.CreateNewCondition("LineCondition3D2N", i+1, [i+1, i+2], prop).SetValue(Kratos.DENSITY, i % 2)

    def test_Union(self):
        mask_1 = Kratos.Expression.NodalExpression(self.model_part)
        mask_2 = Kratos.Expression.NodalExpression(self.model_part)

        Kratos.Expression.VariableExpressionIO.Read(mask_1, Kratos.PRESSURE, False)
        Kratos.Expression.VariableExpressionIO.Read(mask_2, Kratos.DENSITY, False)

        union_mask: Kratos.Expression.NodalExpression = KratosSI.MaskUtils.Union(mask_1, mask_2, 2)
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

        union_mask: Kratos.Expression.NodalExpression = KratosSI.MaskUtils.Intersect(mask_1, mask_2, 2)
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

        union_mask: Kratos.Expression.NodalExpression = KratosSI.MaskUtils.Subtract(mask_1, mask_2, 2)
        Kratos.Expression.VariableExpressionIO.Write(union_mask, Kratos.TEMPERATURE, False)
        for node in self.model_part.Nodes:
            if (node.Id % 3 >= 2 and node.Id % 5 < 2):
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 2)
            else:
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 0)

    def test_GetMaskSize(self):
       mask = Kratos.Expression.NodalExpression(self.model_part)
       Kratos.Expression.VariableExpressionIO.Read(mask, Kratos.DENSITY, False)

       self.assertEqual(KratosSI.MaskUtils.GetMaskSize(mask), self.n - self.n // 5)
       self.assertEqual(KratosSI.MaskUtils.GetMaskSize(mask, 2), self.n - 2 * self.n // 5)
       self.assertEqual(KratosSI.MaskUtils.GetMaskSize(mask, 3), self.n - 3 * self.n // 5)

    def test_GetMaskNoThreshold(self):
        values = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(values, Kratos.HEAT_FLUX, False)

        mask = KratosSI.MaskUtils.GetMask(values)
        Kratos.Expression.VariableExpressionIO.Write(mask, Kratos.TEMPERATURE, False)

        mask_threshold = KratosSI.MaskUtils.GetMaskThreshold(values)
        mask_1 = KratosSI.MaskUtils.GetMask(values, mask_threshold)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(mask - mask_1), 0.0)

        for node in self.model_part.Nodes:
            self.assertEqual(node.GetValue(Kratos.TEMPERATURE), node.Id % 2)

    def test_GetMaskWithhreshold(self):
        values = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(values, Kratos.HEAT_FLUX, False)

        mask = KratosSI.MaskUtils.GetMask(values, 1.0)
        Kratos.Expression.VariableExpressionIO.Write(mask, Kratos.TEMPERATURE, False)

        for node in self.model_part.Nodes:
            self.assertEqual(node.GetValue(Kratos.TEMPERATURE), node.Id % 2)

    def test_Scale(self):
        values = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(values, Kratos.HEAT_FLUX, False)

        float_exp_np = np.arange(0, self.n, dtype=np.float64)
        float_exp = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.CArrayExpressionIO.Read(float_exp, float_exp_np)
        scaled_exp_sum = np.sum(KratosSI.MaskUtils.Scale(float_exp, values).Evaluate())
        self.assertEqual(scaled_exp_sum, (self.n // 2 - 1) * (self.n // 2))

    def test_ClusterMasks(self):
        masks_list: 'list[Kratos.Expression.NodalExpression]' = []
        number_of_masks = 5
        for i in range(number_of_masks):
            mask = np.zeros((self.n), dtype=np.int32)
            for j in range(self.n):
                mask[j] = j % (i + 1)
            mask_exp = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.CArrayExpressionIO.Read(mask_exp, mask)
            masks_list.append(mask_exp)

        cluster_data: 'list[tuple[list[int], ExpressionUnionType]]' = KratosSI.MaskUtils.ClusterMasks(masks_list)
        reference_cluster_indices = [[], [1,2,3,4], [2,3,4], [1,3,4], [2,4], [1,2,3], [3,4]]
        reference_cluster_masks = np.array([
            [1,0,0,0,0,0,0,0,0,0],
            [0,1,0,0,0,0,0,1,0,0],
            [0,0,1,0,0,0,0,0,0,0],
            [0,0,0,1,0,0,0,0,0,1],
            [0,0,0,0,1,0,0,0,1,0],
            [0,0,0,0,0,1,0,0,0,0],
            [0,0,0,0,0,0,1,0,0,0]
        ])

        for i, (cluster_indices, cluster_mask) in enumerate(cluster_data):
            self.assertEqual(cluster_indices, reference_cluster_indices[i])
            self.assertAlmostEqual(np.linalg.norm(cluster_mask.Evaluate() - reference_cluster_masks[i, :]), 0.0, 12)

    def test_GetMasksDividingReferenceMask(self):
        masks_list: 'list[Kratos.Expression.NodalExpression]' = []
        number_of_masks = 5
        for i in range(number_of_masks):
            mask = np.zeros((self.n), dtype=np.int32)
            for j in range(self.n):
                mask[j] = j % (i + 1)
            mask_exp = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.CArrayExpressionIO.Read(mask_exp, mask)
            masks_list.append(mask_exp)

        ref_mask = np.zeros((self.n), dtype=np.int32)
        ref_mask[3:5] = 1
        ref_mask_exp = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.CArrayExpressionIO.Read(ref_mask_exp, ref_mask)
        indices = KratosSI.MaskUtils.GetMasksDividingReferenceMask(ref_mask_exp, masks_list)

        self.assertEqual(indices, [1,2,3])

    def test_FillModelPartUsingClusterMask_Nodes(self):
        exp = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(exp, Kratos.DENSITY, False)
        nodes_mp = self.model_part.CreateSubModelPart("nodal_cluster")
        KratosSI.MaskUtils.FillModelPartUsingClusterMask(nodes_mp, exp)

        for node in self.model_part.Nodes:
            self.assertEqual(nodes_mp.HasNode(node.Id), node.GetValue(Kratos.DENSITY) >= 1.0)

    def test_FillModelPartUsingClusterMask_Elements(self):
        exp = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(exp, Kratos.DENSITY)
        element_mp = self.model_part.CreateSubModelPart("element_cluster")
        KratosSI.MaskUtils.FillModelPartUsingClusterMask(element_mp, exp)

        for element in self.model_part.Elements:
            if element.GetValue(Kratos.DENSITY) >= 1.0:
                self.assertTrue(element_mp.HasElement(element.Id))
                for node in element.GetGeometry():
                    self.assertTrue(element_mp.HasNode(node.Id))
            else:
                self.assertFalse(element_mp.HasElement(element.Id))


if __name__ == '__main__':
    UnitTest.main()
