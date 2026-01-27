
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestEigenvalueSystem(KratosUnittest.TestCase):

    def test_EigenvalueSystemConstruction(self):
        # Create matrices
        G = KratosMultiphysics.SparseContiguousRowGraph(2)
        G.AddEntry(0,0); G.AddEntry(1,1); G.Finalize()

        K = KratosMultiphysics.CsrMatrix(G)
        K.BeginAssemble(); K.AssembleEntry(2.0,0,0); K.AssembleEntry(2.0,1,1); K.FinalizeAssemble()

        M = KratosMultiphysics.CsrMatrix(G)
        M.BeginAssemble(); M.AssembleEntry(1.0,0,0); M.AssembleEntry(1.0,1,1); M.FinalizeAssemble()

        # Vector for eigenvalues
        eig_vals = KratosMultiphysics.SystemVector(2)

        # Dense matrix for eigenvectors
        # KratosMultiphysics.Matrix corresponds to DenseMatrix<double>
        eig_vecs = KratosMultiphysics.Matrix(2, 2)

        # Constructor
        es = KratosMultiphysics.Future.EigenvalueSystem(K, M, eig_vals, eig_vecs, "EigenSystem")
        self.assertEqual(es.Name(), "EigenSystem")

        # Check accessors
        self.assertEqual(es.GetStiffnessMatrix().size1(), 2)
        self.assertEqual(es.GetMassMatrix().size1(), 2)
        self.assertEqual(es.GetEigenvalues().Size(), 2)
        self.assertEqual(es.GetEigenvectors().Size1(), 2)

if __name__ == '__main__':
    KratosUnittest.main()
