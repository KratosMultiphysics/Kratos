import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestSparseMatrixLinearOperator(KratosUnittest.TestCase):
    def test_SparseMatrixLinearOperator(self):
        # Create a simple 2x2 CsrMatrix
        # [ 2.0  3.0 ]
        # [ 1.0  2.0 ]
        G = KratosMultiphysics.SparseContiguousRowGraph(2)
        G.AddEntry(0,0)
        G.AddEntry(0,1)
        G.AddEntry(1,0)
        G.AddEntry(1,1)
        G.Finalize()

        A = KratosMultiphysics.CsrMatrix(G)
        A.BeginAssemble()
        A.AssembleEntry(2.0,0,0)
        A.AssembleEntry(3.0,0,1)
        A.AssembleEntry(1.0,1,0)
        A.AssembleEntry(2.0,1,1)
        A.FinalizeAssemble()

        # Create the CSR matrix linear operator
        lin_op = KratosMultiphysics.Future.SparseMatrixLinearOperator(A)

        # Check operator properties
        self.assertEqual(lin_op.NumRows(), 2)
        self.assertEqual(lin_op.NumCols(), 2)
        self.assertFalse(lin_op.IsMatrixFree())

        # Test GetMatrix
        # Note: GetMatrix returns a reference to the CsrMatrix.
        A_ref = lin_op.GetMatrix()
        self.assertEqual(A_ref.size1(), 2)

        # Set input vector
        x = KratosMultiphysics.SystemVector(2)
        x.SetValue(1.0)

        # Test SpMV
        y = KratosMultiphysics.SystemVector(2)
        y_data = y.Data()
        lin_op.SpMV(x, y)
        self.assertAlmostEqual(y_data[0], 5.0)
        self.assertAlmostEqual(y_data[1], 3.0)

        # Test TransposeSpMV
        y_trans = KratosMultiphysics.SystemVector(2)
        y_trans_data = y_trans.Data()
        lin_op.TransposeSpMV(x, y_trans)
        self.assertAlmostEqual(y_trans_data[0], 3.0)
        self.assertAlmostEqual(y_trans_data[1], 5.0)

        # Test Clear
        lin_op.Clear()
        self.assertEqual(lin_op.NumRows(), 0)
        self.assertEqual(lin_op.NumCols(), 0)

if __name__ == '__main__':
    KratosUnittest.main()
