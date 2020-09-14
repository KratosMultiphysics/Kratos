import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import KratosMultiphysics.python_linear_solver_factory as python_linear_solver_factory
import os

class TestMonotonicityPreservingSolver(KratosUnittest.TestCase):
    ## This case represent
    def test_MonotonicityPreservingSolver(self):
        settings = KratosMultiphysics.Parameters("""
            {
                "solver_type" : "monotonicity_preserving",
                "inner_solver_settings" : {
                    "preconditioner_type": "amg",
                    "solver_type": "amgcl",
                    "smoother_type": "ilu0",
                    "krylov_type": "lgmres",
                    "coarsening_type": "aggregation",
                    "max_iteration": 500,
                    "provide_coordinates": false,
                    "gmres_krylov_space_dimension": 25,
                    "verbosity": 1,
                    "tolerance": 1e-12,
                    "scaling": true,
                    "block_size": 1,
                    "use_block_matrices_if_possible": false,
                    "coarse_enough": 1000,
                    "max_levels": -1,
                    "pre_sweeps": 1,
                    "post_sweeps": 1,
                    "use_gpgpu": false
                }
            }
            """)
        space = KratosMultiphysics.UblasSparseSpace()
        A = KratosMultiphysics.CompressedMatrix()
        this_file_dir = os.path.dirname(os.path.realpath(__file__))
        matrix_file_path = os.path.join(
            this_file_dir,
            "auxiliar_files_for_python_unittest",
            "sparse_matrix_files",
            "A_monotonicity_solver_test.mm")
        rhs_file_path = os.path.join(
            this_file_dir,
            "auxiliar_files_for_python_unittest",
            "sparse_matrix_files",
            "b_monotonicity_solver_test.mm.rhs")
        KratosMultiphysics.ReadMatrixMarketMatrix(matrix_file_path,A)
        b = KratosMultiphysics.Vector()
        KratosMultiphysics.ReadMatrixMarketVector(rhs_file_path,b)
        n = A.Size1()
        x = KratosMultiphysics.Vector(n)
        for dof in x:
            dof = 293.15
        linear_solver = python_linear_solver_factory.ConstructSolver(settings)
        linear_solver.Solve(A,x,b)
        for dof in x:
            self.assertTrue(dof > 0.0)
        solution_norm = space.TwoNorm(x)
        self.assertAlmostEqual(solution_norm, 5515.662307638862)


if __name__ == '__main__':
    KratosUnittest.main()