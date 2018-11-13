from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *

class TestFlags(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Model()
        self.model_part = self.model.CreateModelPart("TestModelPart")
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)

    def tearDown(self):
        pass

    def testFlagSetReset(self):
        node = self.model_part.GetNode(1)

        node.Set(STRUCTURE)
        node.Set(INTERFACE, True)
        node.Set(VISITED, False)

        self.assertTrue(node.IsDefined(STRUCTURE))
        self.assertTrue(node.IsDefined(INTERFACE))
        self.assertTrue(node.IsDefined(VISITED))
        self.assertFalse(node.IsDefined(RIGID))

        self.assertFalse(node.IsNotDefined(STRUCTURE))
        self.assertTrue(node.IsNotDefined(RIGID))

        self.assertTrue(node.Is(STRUCTURE))
        self.assertTrue(node.Is(INTERFACE))
        self.assertFalse(node.Is(VISITED))
        self.assertFalse(node.Is(RIGID))

        self.assertFalse(node.IsNot(STRUCTURE))
        self.assertTrue(node.IsNot(RIGID)) # unset flags default to False


        node.Reset(STRUCTURE)

        self.assertFalse(node.IsDefined(STRUCTURE))

        node.Clear() # reset all flags

        self.assertFalse(node.IsDefined(INTERFACE))
        self.assertFalse(node.IsDefined(VISITED))
        self.assertFalse(node.IsDefined(RIGID))

    def testFlagOr(self):
        node = self.model_part.GetNode(1)

        # set using OR
        node.Set(BOUNDARY | FLUID)

        self.assertTrue(node.Is(BOUNDARY))
        self.assertTrue(node.Is(FLUID))

        # check using OR
        self.assertTrue(node.Is(BOUNDARY))
        self.assertFalse(node.Is(TO_SPLIT))
        self.assertTrue(node.Is(BOUNDARY | TO_SPLIT))

    def testFlagAnd(self):
        node1 = self.model_part.GetNode(1)
        node2 = self.model_part.CreateNewNode(2, 1.0, 1.0, 1.0)

        # the AND of two flags is always nothing
        node1.Set(TO_SPLIT & TO_ERASE)

        self.assertTrue(node1.IsNotDefined(TO_SPLIT))
        self.assertTrue(node1.IsNotDefined(TO_ERASE))

        # Using AND to check for intersection
        node1.Set(ACTIVE | MODIFIED)
        node2.Set(ACTIVE | MPI_BOUNDARY)

        self.assertTrue( (node1 & node2).Is(ACTIVE) )
        self.assertFalse( (node1 & node2).Is(MPI_BOUNDARY) )

    def testFlagFlip(self):
        node = self.model_part.GetNode(1)

        node.Set(RIGID, False)
        node.Flip(RIGID)

        self.assertTrue(node.IsDefined(RIGID))
        self.assertTrue(node.Is(RIGID))

        # unset flags can be flipped too
        self.assertFalse(node.IsDefined(MPI_BOUNDARY))
        node.Flip(MPI_BOUNDARY)

        self.assertTrue(node.IsDefined(MPI_BOUNDARY))
        self.assertTrue(node.Is(MPI_BOUNDARY))

if __name__ == "__main__":
    KratosUnittest.main()