
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestLinearSystem(KratosUnittest.TestCase):
    def test_LinearSystemError(self):
        # Default constructor
        ls = KratosMultiphysics.Future.LinearSystem()

        # Check the empty linear system
        with self.assertRaisesRegex(RuntimeError, "Error: Requested matrix is not initialized for tag 0"):
            ls.GetMatrix(KratosMultiphysics.Future.LinearSystemTag.LHS)

        with self.assertRaisesRegex(RuntimeError, "Error: Requested linear operator is not initialized for tag 0"):
            ls.GetLinearOperator(KratosMultiphysics.Future.LinearSystemTag.LHS)

        with self.assertRaisesRegex(RuntimeError, "Error: Requested vector is not initialized for tag 1"):
            ls.GetVector(KratosMultiphysics.Future.LinearSystemTag.RHS)

        with self.assertRaisesRegex(RuntimeError, "Error: Requested vector is not initialized for tag 2"):
            ls.GetVector(KratosMultiphysics.Future.LinearSystemTag.Dx)

        with self.assertRaisesRegex(RuntimeError, "Error: Requested matrix is not initialized for tag 3"):
            ls.GetMatrix(KratosMultiphysics.Future.LinearSystemTag.M)

        with self.assertRaisesRegex(RuntimeError, "Error: Requested matrix is not initialized for tag 4"):
            ls.GetMatrix(KratosMultiphysics.Future.LinearSystemTag.K)

        with self.assertRaisesRegex(RuntimeError, "Error: Requested matrix is not initialized for tag 5"):
            ls.GetMatrix(KratosMultiphysics.Future.LinearSystemTag.C)

        self.assertEqual(ls.Name, "")

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

        # Check accessors
        self.assertEqual(ls.Name, "TestSystem")
        self.assertEqual(ls.GetMatrix(KratosMultiphysics.Future.LinearSystemTag.LHS).size1(), 2)
        self.assertEqual(ls.GetVector(KratosMultiphysics.Future.LinearSystemTag.RHS).Size(), 2)
        self.assertEqual(ls.GetVector(KratosMultiphysics.Future.LinearSystemTag.Dx).Size(), 2)

        # Check internal LinearOperator creation
        lin_op = ls.GetLinearOperator(KratosMultiphysics.Future.LinearSystemTag.LHS)
        self.assertEqual(lin_op.NumRows, 2)
        self.assertFalse(lin_op.IsMatrixFree)

    def test_LinearSystemWithLinearOperator(self):
        # Create matrix for operator
        G = KratosMultiphysics.SparseContiguousRowGraph(2)
        G.AddEntry(0,0); G.AddEntry(1,1); G.Finalize()
        A = KratosMultiphysics.CsrMatrix(G)
        A.BeginAssemble(); A.AssembleEntry(1.0,0,0); A.AssembleEntry(1.0,1,1); A.FinalizeAssemble()

        # Create vectors
        b = KratosMultiphysics.SystemVector(2)
        b.BeginAssemble()
        b.Assemble([1.0,1.0], [0,1])
        b.FinalizeAssemble()
        x = KratosMultiphysics.SystemVector(2)

        # Create a matrix-based linear operator
        csr_lin_op = KratosMultiphysics.Future.SparseMatrixLinearOperator(A)

        # Constructor
        ls = KratosMultiphysics.Future.LinearSystem(csr_lin_op, b, x, "CsrLinearOperatorTestSystem")
        self.assertEqual(ls.Name, "CsrLinearOperatorTestSystem")
        del csr_lin_op # Manually delete the linear operator to test the ownership transfer

        # Check accessors
        lin_op = ls.GetLinearOperator(KratosMultiphysics.Future.LinearSystemTag.LHS)
        self.assertEqual(lin_op.NumRows, 2)
        self.assertFalse(lin_op.IsMatrixFree)

        # Getting the left hand side matrix should be possible as this is a matrix-based linear operator
        lhs = ls.GetMatrix(KratosMultiphysics.Future.LinearSystemTag.LHS)
        self.assertEqual(lhs.Size1(), 2)
        self.assertEqual(lhs.Size2(), 2)

        # Modify the matrix and apply the corresponding linear operator SpMv product to check that the linear operator is updated accordingly (no internal copy)
        lhs.BeginAssemble()
        lhs.AssembleEntry(2.0,0,0)
        lhs.FinalizeAssemble()
        lin_op.SpMV(b, x)
        self.assertEqual(x[0], 3.0)
        self.assertEqual(x[1], 1.0)

    def test_LinearSystemWithMatrixFreeLinearOperator(self):
        # Create vectors
        b = KratosMultiphysics.SystemVector(2)
        x = KratosMultiphysics.SystemVector(2)

        # Create a matrix-free linear operator
        lin_op = KratosMultiphysics.Future.LinearOperator((2,2))

        # Constructor
        ls = KratosMultiphysics.Future.LinearSystem(lin_op, b, x, "LinearOperatorTestSystem")
        self.assertEqual(ls.Name, "LinearOperatorTestSystem")
        del lin_op # Manually delete the linear operator to test the ownership transfer

        # Check accessors
        lin_op = ls.GetLinearOperator(KratosMultiphysics.Future.LinearSystemTag.LHS)
        self.assertEqual(lin_op.NumRows, 2)
        self.assertTrue(lin_op.IsMatrixFree)

        # GetLeftHandSide should throw as this is pure matrix-free
        with self.assertRaisesRegex(RuntimeError, "Linear operator referring to matrix 0 is matrix free"):
            ls.GetMatrix(KratosMultiphysics.Future.LinearSystemTag.LHS)

if __name__ == '__main__':
    KratosUnittest.main()
