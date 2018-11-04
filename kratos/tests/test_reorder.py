from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

import os
import sys


def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class TestReorder(KratosUnittest.TestCase):

    def setUp(self):
        if (sys.version_info < (3, 2)):
            self.assertRaisesRegex = self.assertRaisesRegexp

    def test_reorder(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.CreateNewNode(1,0,0,0)
        model_part.CreateNewNode(2,1,0,0)
        model_part.CreateNewNode(3,2,0,0)
        model_part.CreateNewNode(4,3,0,0)
        model_part.CreateNewNode(5,0,1,0)
        model_part.CreateNewNode(6,1,1,0)
        model_part.CreateNewNode(7,2,1,0)
        model_part.CreateNewNode(8,3,1,0)
        model_part.CreateNewNode(9,0,2,0)
        model_part.CreateNewNode(10,1,2,0)
        model_part.CreateNewNode(11,2,2,0)
        model_part.CreateNewNode(12,3,2,0)
        #nodes are numbered as
        #9 10 11 12
        #5 6  7  8
        #1 2  3  4
        model_part.CreateNewElement("Element2D4N",1,[1,2,6,5], model_part.GetProperties()[1])
        model_part.CreateNewElement("Element2D4N",2,[2,3,7,6], model_part.GetProperties()[1])
        model_part.CreateNewElement("Element2D4N",3,[3,4,8,7], model_part.GetProperties()[1])
        model_part.CreateNewElement("Element2D4N",4,[5,6,10,9], model_part.GetProperties()[1])
        model_part.CreateNewElement("Element2D4N",5,[6,7,11,10], model_part.GetProperties()[1])
        model_part.CreateNewElement("Element2D4N",6,[7,8,12,11], model_part.GetProperties()[1])

        tmp = KratosMultiphysics.Parameters("{}")
        KratosMultiphysics.ReorderAndOptimizeModelPartProcess(model_part,tmp).Execute()

        #output is
        #7   8 11 12
        #4   5  6 10
        #1   2  3 9
        self.assertEqual(model_part.Nodes[1].X, 0.0);  self.assertEqual(model_part.Nodes[1].Y, 0.0)
        self.assertEqual(model_part.Nodes[2].X, 1.0);  self.assertEqual(model_part.Nodes[2].Y, 0.0)
        self.assertEqual(model_part.Nodes[3].X, 2.0);  self.assertEqual(model_part.Nodes[3].Y, 0.0)
        self.assertEqual(model_part.Nodes[9].X, 3.0);  self.assertEqual(model_part.Nodes[9].Y, 0.0)
        self.assertEqual(model_part.Nodes[4].X, 0.0);  self.assertEqual(model_part.Nodes[4].Y, 1.0)
        self.assertEqual(model_part.Nodes[5].X, 1.0);  self.assertEqual(model_part.Nodes[5].Y, 1.0)
        self.assertEqual(model_part.Nodes[6].X, 2.0);  self.assertEqual(model_part.Nodes[6].Y, 1.0)
        self.assertEqual(model_part.Nodes[10].X, 3.0);  self.assertEqual(model_part.Nodes[10].Y, 1.0)
        self.assertEqual(model_part.Nodes[7].X, 0.0);  self.assertEqual(model_part.Nodes[7].Y, 2.0)
        self.assertEqual(model_part.Nodes[8].X, 1.0);  self.assertEqual(model_part.Nodes[8].Y, 2.0)
        self.assertEqual(model_part.Nodes[11].X, 2.0);  self.assertEqual(model_part.Nodes[11].Y, 2.0)
        self.assertEqual(model_part.Nodes[12].X, 3.0);  self.assertEqual(model_part.Nodes[12].Y, 2.0)

        #for node in model_part.Nodes:
            #print(node.Id, node.X,node.Y)

        #for elem in model_part.Elements:
            #tmp = []
            #for node in elem.GetNodes():
                #tmp.append(node.Id)
            #print(elem.Id,tmp)


if __name__ == '__main__':
    KratosUnittest.main()
