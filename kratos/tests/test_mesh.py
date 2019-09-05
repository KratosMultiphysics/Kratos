from __future__ import print_function, absolute_import, division

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestMesh(KratosUnittest.TestCase):
    def setUp(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("test")

        for i in range(5):
            self.model_part.CreateNewNode(i+1, i*1.0, 0.0, 0.0)

        self.mesh = self.model_part.GetMesh()

    def test_nodes_interface(self):
        self.assertEqual(self.mesh.NumberOfNodes(), 5)
        self.assertEqual(len(self.mesh.Nodes), 5)

        self.mesh.Nodes.append(self.mesh.Nodes[1])
        self.assertEqual(len(self.mesh.Nodes), 6)

if __name__ == '__main__':
    KratosUnittest.main()
