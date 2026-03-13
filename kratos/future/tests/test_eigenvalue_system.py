
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.Future as KratosFuture

class TestEigenvalueSystem(KratosUnittest.TestCase):
    def test_eigenvalue_system_construction(self):
        # Default constructor
        es = KratosFuture.EigenvalueSystem()
        self.assertEqual(es.Name(), "")

    def test_eigenvalue_system_full_construction(self):
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
        # Note: signature expects DenseMatrixPointerType which is shared_ptr<DenseMatrix>.
        # Python binding for Matrix usually implies value/ref logic, but let's see if it accepts the object.
        # If it requires a pointer/shared_ptr explicitly wrapper, it might fail if not handled by pybind11 automatically.
        # However, for Matrix, Kratos usually doesn't use shared_ptr in python interface often, but here C++ demands it.
        # If this fails, we might need a workaround or check how Matrix pointers are handled.
        es = KratosFuture.EigenvalueSystem(K, M, eig_vals, eig_vecs, "EigenSystem")

        self.assertEqual(es.Name(), "EigenSystem")

        self.assertEqual(es.GetStiffnessMatrix().size1(), 2)
        self.assertEqual(es.GetMassMatrix().size1(), 2)
        self.assertEqual(es.GetEigenvalues().Size(), 2)
        self.assertEqual(es.GetEigenvectors().Size1(), 2)

if __name__ == '__main__':
    KratosUnittest.main()
