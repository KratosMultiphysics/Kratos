
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestLinearSystem(KratosUnittest.TestCase):
    def test_LinearSystemError(self):
        # Default constructor
        ls = KratosMultiphysics.Future.LinearSystem()
        # IsMatrixFree throws if not initialized
        with self.assertRaisesRegex(RuntimeError, "Linear system has no matrix or linear operator"):
            ls.IsMatrixFree()

        self.assertEqual(ls.Name(), "")

    def test_LinearSystemWithMatrix(self):
        # Create matrix
        G = KratosMultiphysics.SparseContiguousRowGraph(2)
        G.AddEntry(0,0); G.AddEntry(1,1); G.Finalize()
        A = KratosMultiphysics.CsrMatrix(G)
        A.BeginAssemble(); A.AssembleEntry(1.0,0,0); A.AssembleEntry(1.0,1,1); A.FinalizeAssemble()

        # Create vectors
        b = KratosMultiphysics.SystemVector(2)
        x = KratosMultiphysics.SystemVector(2)

        # Constructor
        ls = KratosMultiphysics.Future.LinearSystem(A, b, x, "TestSystem")
        self.assertEqual(ls.Name(), "TestSystem")
        self.assertFalse(ls.IsMatrixFree())

        # Check accessors
        self.assertEqual(ls.GetLeftHandSide().size1(), 2)
        self.assertEqual(ls.GetRightHandSide().Size(), 2)
        self.assertEqual(ls.GetSolution().Size(), 2)

        # Check internal LinearOperator creation
        lin_op = ls.GetLinearOperator()
        self.assertEqual(lin_op.NumRows(), 2)
        self.assertFalse(lin_op.IsMatrixFree())

    def test_LinearSystemWithLinearOperator(self):
        # Create matrix for operator
        G = KratosMultiphysics.SparseContiguousRowGraph(2)
        G.AddEntry(0,0); G.AddEntry(1,1); G.Finalize()
        A = KratosMultiphysics.CsrMatrix(G)
        A.BeginAssemble(); A.AssembleEntry(1.0,0,0); A.AssembleEntry(1.0,1,1); A.FinalizeAssemble()

        lin_op = KratosMultiphysics.Future.SparseMatrixLinearOperator(A)

        # Create vectors
        b = KratosMultiphysics.SystemVector(2)
        x = KratosMultiphysics.SystemVector(2)

        # Constructor
        ls = KratosMultiphysics.Future.LinearSystem(lin_op, b, x, "LinearOperatorTestSystem")
        self.assertEqual(ls.Name(), "LinearOperatorTestSystem")
        self.assertTrue(ls.IsMatrixFree())

        # Check accessors
        self.assertEqual(ls.GetLinearOperator().NumRows(), 2)

        # GetLeftHandSide should throw if no matrix
        with self.assertRaisesRegex(RuntimeError, "Left-hand side matrix is not initialized"):
            ls.GetLeftHandSide()

if __name__ == '__main__':
    KratosUnittest.main()
