from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *

class TestFlags(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Model()
        self.model_part = self.model.CreateModelPart("TestModelPart")
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)

    def test_CopyConstructor(self):
        flag1 = ACTIVE | STRUCTURE | MPI_BOUNDARY
        flag2 = Flags(flag1)

        flag1 |= VISITED | RIGID # modify the original flags, this should not affect the copy-constructed flags
        self.assertTrue(flag1.Is(VISITED))
        self.assertTrue(flag1.Is(RIGID))

        self.assertTrue(flag2.Is(ACTIVE))
        self.assertTrue(flag2.Is(STRUCTURE))
        self.assertTrue(flag2.Is(MPI_BOUNDARY))
        self.assertFalse(flag2.Is(VISITED))
        self.assertFalse(flag2.Is(RIGID))

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

        # the AND of two named flags sets both
        node1.Set(TO_SPLIT & TO_ERASE)

        # but they are set to false, because the bit on the other argument is unset (undefined) & true = (undefined)
        self.assertFalse(node1.Is(TO_SPLIT))
        self.assertFalse(node1.Is(TO_ERASE))

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

    def testFlagAndEqual(self):
        #       both true | both false | opposite sets | first true   | first false | second true | second false
        flag1 = ACTIVE    | (RIGID).AsFalse()  | STRUCTURE     | MPI_BOUNDARY | (PERIODIC).AsFalse()
        flag2 = ACTIVE    | (RIGID).AsFalse()  | (STRUCTURE).AsFalse() |                              INLET       | (OUTLET).AsFalse()

        flag1 &= flag2
        # true (defined) & true (defined) = true (defined)
        self.assertTrue(flag1.IsDefined(ACTIVE))
        self.assertTrue(flag1.Is(ACTIVE))
        # false (defined) & false (defined) = false (defined)
        self.assertTrue(flag1.IsDefined(RIGID))
        self.assertFalse(flag1.Is(RIGID))
        # true (defined) & false (defined) = false (defined)
        self.assertTrue(flag1.IsDefined(STRUCTURE))
        self.assertFalse(flag1.Is(STRUCTURE))
        # true (defined) & (undefined) = false (defined) (undefined defaults to false)
        self.assertTrue(flag1.IsDefined(MPI_BOUNDARY))
        self.assertFalse(flag1.Is(MPI_BOUNDARY))
        # false (defined) & (undefined) = false (defined) (undefined defaults to false)
        self.assertTrue(flag1.IsDefined(PERIODIC))
        self.assertFalse(flag1.Is(PERIODIC))
        # (undefined) & true (defined) = false (defined) (undefined defaults to false)
        self.assertTrue(flag1.IsDefined(INLET))
        self.assertFalse(flag1.Is(INLET))
        # (undefied) & false (defined) = false (defined) (undefined defaults to false)
        self.assertTrue(flag1.IsDefined(OUTLET))
        self.assertFalse(flag1.Is(OUTLET))
        # (undefined) & (undefined) = (undefined)
        self.assertFalse(flag1.IsDefined(ISOLATED))
        self.assertFalse(flag1.Is(ISOLATED))

    def testFlagOrEqual(self):
        #       both true | both false | opposite sets | first true   | first false | second true | second false
        flag1 = ACTIVE    | (RIGID).AsFalse()  | STRUCTURE     | MPI_BOUNDARY | (PERIODIC).AsFalse()
        flag2 = ACTIVE    | (RIGID).AsFalse()  | (STRUCTURE).AsFalse() |                              INLET       | (OUTLET).AsFalse()

        flag1 |= flag2
        # true (defined) | true (defined) = true (defined)
        self.assertTrue(flag1.IsDefined(ACTIVE))
        self.assertTrue(flag1.Is(ACTIVE))
        # false (defined) | false (defined) = false (defined)
        self.assertTrue(flag1.IsDefined(RIGID))
        self.assertFalse(flag1.Is(RIGID))
        # true (defined) | false (defined) = true (defined)
        self.assertTrue(flag1.IsDefined(STRUCTURE))
        self.assertTrue(flag1.Is(STRUCTURE))
        # true (defined) | (undefined) = true (defined) (undefined defaults to false)
        self.assertTrue(flag1.IsDefined(MPI_BOUNDARY))
        self.assertTrue(flag1.Is(MPI_BOUNDARY))
        # false (defined) | (undefined) = false (defined) (undefined defaults to false)
        self.assertTrue(flag1.IsDefined(PERIODIC))
        self.assertFalse(flag1.Is(PERIODIC))
        # (undefined) | true (defined) = true (defined) (undefined defaults to false)
        self.assertTrue(flag1.IsDefined(INLET))
        self.assertTrue(flag1.Is(INLET))
        # (undefied) | false (defined) = false (defined) (undefined defaults to false)
        self.assertTrue(flag1.IsDefined(OUTLET))
        self.assertFalse(flag1.Is(OUTLET))
        # (undefined) | (undefined) = (undefined)
        self.assertFalse(flag1.IsDefined(ISOLATED))
        self.assertFalse(flag1.Is(ISOLATED))

if __name__ == "__main__":
    KratosUnittest.main()
