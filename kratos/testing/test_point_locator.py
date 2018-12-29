from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

class TestPointLocator(KratosUnittest.TestCase):

    def test_point_locator_2d(self):
        current_model = KratosMultiphysics.Model()
        model_part1 = current_model.CreateModelPart("Main")

        # 3 - 4
        # | / |
        # 1 - 2
        model_part1.CreateNewNode(1,0.0,0.0,0.0)
        model_part1.CreateNewNode(2,1.0,0.0,0.0)
        model_part1.CreateNewNode(3,0.0,1.0,0.0)
        model_part1.CreateNewNode(4,1.0,1.0,0.0)

        model_part1.CreateNewElement("Element2D3N", 1, [1,2,4], model_part1.GetProperties()[1])
        model_part1.CreateNewElement("Element2D3N", 2, [1,4,3], model_part1.GetProperties()[1])

        locator = KratosMultiphysics.BinBasedFastPointLocator2D(model_part1)
        pos = KratosMultiphysics.Array3()
        pos[0] = 2.0/3.0
        pos[1] = 1.0/3.0
        pos[2] = 0.0

        [found,N,pelem] = locator.FindPointOnMesh(pos,1e-9,1000)

        self.assertTrue(found)
        self.assertEqual(N[0], 1.0/3.0)
        self.assertEqual(N[1], 1.0/3.0)
        self.assertEqual(N[2], 1.0/3.0)
        self.assertEqual(pelem.Id, 1)

    def test_point_locator_3d(self):
        current_model = KratosMultiphysics.Model()
        model_part1 = current_model.CreateModelPart("Main")

        model_part1.CreateNewNode(1,0.0,0.0,0.0)
        model_part1.CreateNewNode(2,1.0,0.0,0.0)
        model_part1.CreateNewNode(3,0.0,1.0,0.0)
        model_part1.CreateNewNode(4,0.0,0.0,1.0)

        model_part1.CreateNewElement("Element3D4N", 1, [1,2,3,4], model_part1.GetProperties()[1])

        locator = KratosMultiphysics.BinBasedFastPointLocator3D(model_part1)
        pos = KratosMultiphysics.Array3()
        pos[0] = 0.25
        pos[1] = 0.25
        pos[2] = 0.25

        [found,N,pelem] = locator.FindPointOnMesh(pos,1e-9,1000)

        self.assertTrue(found)
        self.assertEqual(N[0], 0.25)
        self.assertEqual(N[1], 0.25)
        self.assertEqual(N[2], 0.25)
        self.assertEqual(N[3], 0.25)
        self.assertEqual(pelem.Id, 1)

if __name__ == '__main__':
    KratosUnittest.main()
