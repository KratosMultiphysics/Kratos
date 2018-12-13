from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
from vtk_output_process import VtkOutputProcessPython

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestVtkOutputProcess(KratosUnittest.TestCase):
    def __SetupModelPart(self, model):
        mp = model.CreateModelPart("Main")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        #create nodes
        self.mp.CreateNewNode(1, 0.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(2, 0.00000, 0.50000, 0.00000)
        self.mp.CreateNewNode(3, 0.50000, 1.00000, 0.00000)
        self.mp.CreateNewNode(4, 0.50000, 0.50000, 0.00000)
        self.mp.CreateNewNode(5, 0.00000, 0.00000, 0.00000)
        self.mp.CreateNewNode(6, 1.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(7, 1.00000, 0.50000, 0.00000)
        self.mp.CreateNewNode(8, 0.50000, 0.00000, 0.00000)
        self.mp.CreateNewNode(9, 1.00000, 0.00000, 0.00000)
        self.mp.CreateNewNode(10, 1.50000, 1.00000, 0.00000)
        self.mp.CreateNewNode(11, 1.50000, 0.50000, 0.00000)
        self.mp.CreateNewNode(12, 1.50000, 0.00000, 0.00000)
        self.mp.CreateNewNode(13, 2.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(14, 2.00000, 0.50000, 0.00000)
        self.mp.CreateNewNode(15, 2.00000, 0.00000, 0.00000)
        self.mp.CreateNewNode(16, 1.00000, 1.00000, 0.00000)
        self.mp.CreateNewNode(17, 1.00000, 0.50000, 0.00000)
        self.mp.CreateNewNode(18, 1.00000, 0.00000, 0.00000)

        #create a submodelpart for boundary conditions
        bcs = self.mp.CreateSubModelPart("FixedEdgeNodes")
        bcs.AddNodes([1, 2, 5])

        bcmn = self.mp.CreateSubModelPart("MovingNodes")
        bcmn.AddNodes([13, 14, 15])

        #create Element
        self.mp.CreateNewElement("Element2D4N", 1,
                            [14, 11, 12, 15], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 2,
                            [13, 10, 11, 14], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 3,
                            [11, 17, 18, 12], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 4,
                            [10, 16, 17, 11], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 5,
                            [2, 4, 3, 1], self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 6, [5, 8, 4, 2],
                            self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 7, [4, 7, 6, 3],
                            self.mp.GetProperties()[1])
        self.mp.CreateNewElement("Element2D4N", 8, [8, 9, 7, 4],
                            self.mp.GetProperties()[1])

    def test_vtk_io(self):
        current_model = KratosMultiphysics.Model()