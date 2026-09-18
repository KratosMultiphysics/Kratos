import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

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

    def test_Union(self):
        mask_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        mask_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.DENSITY)

        mask_1.CollectData()
        mask_2.CollectData()

        union_mask = KratosSI.MaskUtils.Union(mask_1, mask_2, 2)
        Kratos.TensorAdaptors.VariableTensorAdaptor(union_mask, Kratos.TEMPERATURE, copy=False).StoreData()
        for node in self.model_part.Nodes:
            if (node.Id % 3 >= 2 or node.Id % 5 >= 2):
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 2)
            else:
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 0)

    def test_Intersect(self):
        mask_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        mask_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.DENSITY)

        mask_1.CollectData()
        mask_2.CollectData()

        intersect_mask = KratosSI.MaskUtils.Intersect(mask_1, mask_2, 2)
        Kratos.TensorAdaptors.VariableTensorAdaptor(intersect_mask, Kratos.TEMPERATURE, copy=False).StoreData()
        for node in self.model_part.Nodes:
            if (node.Id % 3 >= 2 and node.Id % 5 >= 2):
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 2)
            else:
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 0)

    def test_Substract(self):
        mask_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        mask_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.DENSITY)

        mask_1.CollectData()
        mask_2.CollectData()

        subtract_mask = KratosSI.MaskUtils.Subtract(mask_1, mask_2, 2)
        Kratos.TensorAdaptors.VariableTensorAdaptor(subtract_mask, Kratos.TEMPERATURE, copy=False).StoreData()
        for node in self.model_part.Nodes:
            if (node.Id % 3 >= 2 and node.Id % 5 < 2):
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 2)
            else:
                self.assertEqual(node.GetValue(Kratos.TEMPERATURE), 0)

    def test_GetMaskSize(self):
       mask = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.DENSITY)
       mask.CollectData()

       self.assertEqual(KratosSI.MaskUtils.GetMaskSize(mask), self.n - self.n // 5)
       self.assertEqual(KratosSI.MaskUtils.GetMaskSize(mask, 2), self.n - 2 * self.n // 5)
       self.assertEqual(KratosSI.MaskUtils.GetMaskSize(mask, 3), self.n - 3 * self.n // 5)

    def test_GetMaskNoThreshold(self):
        values = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.HEAT_FLUX)
        values.CollectData()

        mask = KratosSI.MaskUtils.GetMask(values)
        Kratos.TensorAdaptors.VariableTensorAdaptor(mask, Kratos.TEMPERATURE, copy=False).StoreData()

        mask_threshold = KratosSI.MaskUtils.GetMaskThreshold(values)
        mask_1 = KratosSI.MaskUtils.GetMask(values, mask_threshold)
        self.assertAlmostEqual(np.linalg.norm(mask.data - mask_1.data), 0.0)

        for node in self.model_part.Nodes:
            self.assertEqual(node.GetValue(Kratos.TEMPERATURE), node.Id % 2)

    def test_GetMaskWithThreshold(self):
        values = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.HEAT_FLUX)
        values.CollectData()

        mask = KratosSI.MaskUtils.GetMask(values, 1.0)
        Kratos.TensorAdaptors.VariableTensorAdaptor(mask, Kratos.TEMPERATURE, copy=False).StoreData()

        for node in self.model_part.Nodes:
            self.assertEqual(node.GetValue(Kratos.TEMPERATURE), node.Id % 2)

    def test_Scale(self):
        values = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.HEAT_FLUX)
        values.CollectData()

        float_ta = values.Clone()
        float_ta.data[:] = np.arange(0, self.n, dtype=np.float64)
        scaled_exp_sum = np.sum(KratosSI.MaskUtils.Scale(float_ta, values).data)
        self.assertEqual(scaled_exp_sum, (self.n // 2 - 1) * (self.n // 2))

    def test_ClusterMasks(self):
        masks_list: 'list[Kratos.TensorAdaptors.DoubleTensorAdaptor]' = []
        number_of_masks = 5
        for i in range(number_of_masks):
            mask_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
            for j in range(self.n):
                mask_ta.data[j] = j % (i + 1)
            masks_list.append(mask_ta)

        cluster_data = KratosSI.MaskUtils.ClusterMasks(masks_list)
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
            self.assertAlmostEqual(np.linalg.norm(cluster_mask.data - reference_cluster_masks[i, :]), 0.0, 12)

    def test_GetMasksDividingReferenceMask(self):
        masks_list: 'list[Kratos.TensorAdaptors.DoubleTensorAdaptor]' = []
        number_of_masks = 5
        for i in range(number_of_masks):
            mask_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
            for j in range(self.n):
                mask_ta.data[j] = j % (i + 1)
            masks_list.append(mask_ta)

        ref_mask = np.zeros((self.n), dtype=np.int32)
        ref_mask[3:5] = 1
        ref_mask_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        ref_mask_ta.data[:] = ref_mask
        indices = KratosSI.MaskUtils.GetMasksDividingReferenceMask(ref_mask_ta, masks_list)

        self.assertEqual(indices, [1,2,3])

if __name__ == '__main__':
    UnitTest.main()
