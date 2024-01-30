
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

class TestPointLocator(KratosUnittest.TestCase):

    def test_point_locator_2d(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("Main")

        # 3 - 4
        # | / |
        # 1 - 2
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,0.0,0.0)
        mp.CreateNewNode(3,0.0,1.0,0.0)
        mp.CreateNewNode(4,1.0,1.0,0.0)

        mp.CreateNewElement("Element2D3N", 1, [1,4,3], mp.GetProperties()[1])
        mp.CreateNewElement("Element2D3N", 2, [1,2,4], mp.GetProperties()[1])

        locator = KratosMultiphysics.BinBasedFastPointLocator2D(mp)
        locator.UpdateSearchDatabase()
        pos = KratosMultiphysics.Array3()
        pos[0] = 2.0/3.0
        pos[1] = 1.0/3.0
        pos[2] = 0.0

        [found,N,pelem] = locator.FindPointOnMesh(pos)

        self.assertTrue(found)
        self.assertAlmostEqual(N[0], 1.0/3.0)
        self.assertAlmostEqual(N[1], 1.0/3.0)
        self.assertAlmostEqual(N[2], 1.0/3.0)
        self.assertEqual(pelem.Id, 2)

    def test_point_locator_3d(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("Main")

        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,0.0,0.0)
        mp.CreateNewNode(3,0.0,1.0,0.0)
        mp.CreateNewNode(4,0.0,0.0,1.0)

        mp.CreateNewElement("Element3D4N", 1, [1,2,3,4], mp.GetProperties()[1])

        locator = KratosMultiphysics.BinBasedFastPointLocator3D(mp)
        locator.UpdateSearchDatabase()
        pos = KratosMultiphysics.Array3()
        pos[0] = 0.25
        pos[1] = 0.25
        pos[2] = 0.25

        [found,N,pelem] = locator.FindPointOnMesh(pos)

        self.assertTrue(found)
        self.assertAlmostEqual(N[0], 0.25)
        self.assertAlmostEqual(N[1], 0.25)
        self.assertAlmostEqual(N[2], 0.25)
        self.assertAlmostEqual(N[3], 0.25)
        self.assertEqual(pelem.Id, 1)

if __name__ == '__main__':
    KratosUnittest.main()
