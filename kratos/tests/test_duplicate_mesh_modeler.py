import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestDuplicateMeshModeler(KratosUnittest.TestCase):
    def test_DuplicateMeshModeler(self):
        model = KratosMultiphysics.Model()
        mp1 = model.CreateModelPart("part_1")
        prop0 = mp1.CreateNewProperties(0)
        mp1.CreateNewNode(1, 0, 0 ,0)
        mp1.CreateNewNode(2, 0, 0 ,0)
        mp1.CreateNewNode(3, 0, 0 ,0)
        mp1.CreateNewNode(4, 0, 0 ,0)
        mp1.CreateNewNode(5, 0, 0 ,0)
        mp1.CreateNewElement("Element3D4N", 1, [1,2,3,4], prop0)
        mp1.CreateNewElement("Element3D4N", 2, [1,2,3,5], prop0)
        mp1.CreateNewCondition("SurfaceCondition3D3N", 1, [1,2,4], prop0)
        mp1.CreateNewCondition("SurfaceCondition3D3N", 2, [1,2,5], prop0)
        mp1_1 = mp1.CreateSubModelPart("sub_1")
        mp1_1.AddNodes([1,2,3,4])
        mp1_1.AddElements([1])
        mp1_1.AddConditions([1])

        mp2 = model.CreateModelPart("part_2")
        KratosMultiphysics.DuplicateMeshModeler(mp1).GenerateMesh(
            mp2, "Element3D4N", "SurfaceCondition3D3N")

        self.assertEqual(mp2.NumberOfNodes(), 5)
        self.assertEqual(mp2.NumberOfElements(), 2)
        self.assertEqual(mp2.NumberOfConditions(), 2)
        self.assertTrue(mp2.HasSubModelPart("sub_1"))
        mp_2_1 = mp2.GetSubModelPart("sub_1")
        self.assertEqual(mp_2_1.NumberOfNodes(), 4)
        self.assertEqual(mp_2_1.NumberOfElements(), 1)
        self.assertEqual(mp_2_1.NumberOfConditions(), 1)

if __name__ == "__main__":
    KratosUnittest.main()
