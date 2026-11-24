
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestEigenSparseLUSolver(KratosUnittest.TestCase):
    # Same 5x5 system as the C++ tests in
    # kratos/future/tests/cpp_tests/test_future_linear_solvers.cpp
    def _build_system(self):
        G = KratosMultiphysics.SparseContiguousRowGraph(5)
        for row, col in [
            (0, 0), (0, 3),
            (1, 1), (1, 2),
            (2, 2), (2, 1), (2, 4),
            (3, 3), (3, 0), (3, 4),
            (4, 4), (4, 2), (4, 3)
        ]:
            G.AddEntry(row, col)
        G.Finalize()
        A = KratosMultiphysics.CsrMatrix(G)
        A.BeginAssemble()
        A.AssembleEntry(10.0, 0, 0); A.AssembleEntry(-2.0, 0, 3)
        A.AssembleEntry(8.0, 1, 1);  A.AssembleEntry(-1.0, 1, 2)
        A.AssembleEntry(5.0, 2, 2);  A.AssembleEntry(-1.0, 2, 1); A.AssembleEntry(-2.0, 2, 4)
        A.AssembleEntry(7.0, 3, 3);  A.AssembleEntry(-2.0, 3, 0); A.AssembleEntry(-1.0, 3, 4)
        A.AssembleEntry(6.0, 4, 4);  A.AssembleEntry(-2.0, 4, 2); A.AssembleEntry(-1.0, 4, 3)
        A.FinalizeAssemble()
        b = KratosMultiphysics.SystemVector(5)
        b.BeginAssemble()
        b.Assemble([3.0, 7.0, 2.0, 5.0, 1.0], [0, 1, 2, 3, 4])
        b.FinalizeAssemble()
        return A, b

    def test_parameters_constructor(self):
        # The solver can be constructed from a Parameters object
        solver = KratosMultiphysics.Future.EigenSparseLUSolver(
            KratosMultiphysics.Parameters("""{"solver_type": "eigen_sparse_lu"}"""))
        self.assertIsNotNone(solver)

    def test_solve(self):
        A, b = self._build_system()
        x = KratosMultiphysics.SystemVector(5)
        ls = KratosMultiphysics.Future.LinearSystem(A, b, x, "eigen_lu_system")
        solver = KratosMultiphysics.Future.EigenSparseLUSolver()
        solver.Initialize(ls)
        solver.InitializeSolutionStep(ls)
        self.assertTrue(solver.PerformSolutionStep(ls))
        reference = [0.487946221604, 0.979601298099, 0.836810384794, 0.93973110802, 0.602225312935]
        for i in range(5):
            self.assertAlmostEqual(x[i], reference[i], delta=1.0e-10)

    def test_singular_matrix_raises(self):
        G = KratosMultiphysics.SparseContiguousRowGraph(2)
        G.AddEntry(0, 0); G.AddEntry(0, 1); G.AddEntry(1, 0); G.AddEntry(1, 1)
        G.Finalize()
        A = KratosMultiphysics.CsrMatrix(G)
        A.BeginAssemble()
        A.AssembleEntry(1.0, 0, 0); A.AssembleEntry(1.0, 0, 1)
        A.AssembleEntry(1.0, 1, 0); A.AssembleEntry(1.0, 1, 1)
        A.FinalizeAssemble()
        b = KratosMultiphysics.SystemVector(2)
        b.BeginAssemble()
        b.Assemble([1.0, 1.0], [0, 1])
        b.FinalizeAssemble()
        x = KratosMultiphysics.SystemVector(2)
        ls = KratosMultiphysics.Future.LinearSystem(A, b, x, "singular_system")
        solver = KratosMultiphysics.Future.EigenSparseLUSolver()
        with self.assertRaisesRegex(RuntimeError, "the LHS matrix is singular"):
            solver.InitializeSolutionStep(ls)

if __name__ == '__main__':
    KratosUnittest.main()
