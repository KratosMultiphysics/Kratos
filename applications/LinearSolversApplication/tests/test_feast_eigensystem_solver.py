import os
import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import eigen_solver_factory
import KratosMultiphysics.LinearSolversApplication as LinearSolversApplication

from math import sqrt

@KratosUnittest.skipUnless(LinearSolversApplication.HasFEAST(),"FEAST not found in LinearSolversApplication, skipping.")
class TestFeastEigensystemSolver(KratosUnittest.TestCase):
    def test_real_symmetric_gev(self):

        space = KratosMultiphysics.UblasSparseSpace()

        settings = KratosMultiphysics.Parameters('''{
            "solver_type": "feast",
            "symmetric": true,
            "number_of_eigenvalues": 3,
            "search_lowest_eigenvalues": true,
            "sort_eigenvalues": true,
            "sort_order": "sr",
            "e_min": 0.0,
            "e_max": 0.2,
            "echo_level": 0
         }''')

        K = KratosMultiphysics.CompressedMatrix()

        this_file_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(this_file_dir)))
        matrix_file_path = os.path.join(base_dir, "kratos", "tests", "test_files", "sparse_matrix_files", "A.mm")

        file_read = KratosMultiphysics.ReadMatrixMarketMatrix(matrix_file_path, K) # symmetric test matrix
        self.assertTrue(file_read, msg="The matrix file could not be read")

        n = K.Size1()
        self.assertEqual(n, 900)

        # create an identity matrix
        M = KratosMultiphysics.CompressedMatrix(n, n)
        for i in range(n):
            M[i, i] = 1.0

        # create result containers (they will be resized inside the solver)
        eigenvalues = KratosMultiphysics.Vector(n)
        eigenvectors = KratosMultiphysics.Matrix(n, 1)

        # Construct the solver
        eigen_solver = eigen_solver_factory.ConstructSolver(settings)

        # Solve
        eigen_solver.Solve(K, M, eigenvalues, eigenvectors)

        self.assertAlmostEqual(eigenvalues[0], 0.061463, 4)
        self.assertAlmostEqual(eigenvalues[1], 0.15318431112733275, 7)
        self.assertAlmostEqual(eigenvalues[2], 0.153184311127333, 7)

        for i in range(eigenvalues.Size()):
            eigenvector = KratosMultiphysics.Vector(n)
            for j in range(n):
                eigenvector[j] = eigenvectors[i,j]

            _aux_1 = KratosMultiphysics.Vector(n,0)
            _aux_2 = KratosMultiphysics.Vector(n,0)

            space.Mult(K,eigenvector,_aux_1)
            space.Mult(M,eigenvector,_aux_2)
            _aux_2 *= eigenvalues[i]

            self.assertVectorAlmostEqual(_aux_1, _aux_2)

    def test_real_general_gev(self):

        space = KratosMultiphysics.UblasComplexSparseSpace()

        settings = KratosMultiphysics.Parameters('''{
            "solver_type": "feast",
            "symmetric": false,
            "number_of_eigenvalues": 3,
            "sort_eigenvalues": true,
            "sort_order": "sr",
            "e_mid_re": 10.0,
            "e_mid_im": 0.0,
            "e_r": 3.0,
            "echo_level": 0
         }''')

        K = KratosMultiphysics.CompressedMatrix()

        this_file_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(this_file_dir)))
        matrix_file_path = os.path.join(base_dir, "kratos", "tests", "test_files", "sparse_matrix_files", "small_A.mm")

        file_read = KratosMultiphysics.ReadMatrixMarketMatrix(matrix_file_path, K) # symmetric test matrix
        self.assertTrue(file_read, msg="The matrix file could not be read")

        n = K.Size1()

        # create an identity matrix
        M = KratosMultiphysics.CompressedMatrix(n, n)
        for i in range(n):
            M[i, i] = 1

        # create result containers (they will be resized inside the solver)
        # eigenvalues and vectors of unsymmetric matrices are required to be real here
        eigenvalues = KratosMultiphysics.ComplexVector(n)
        eigenvectors = KratosMultiphysics.ComplexMatrix(n, 1)

        # Construct the solver
        eigen_solver = eigen_solver_factory.ConstructSolver(settings)

        # Solve
        eigen_solver.Solve(K, M, eigenvalues, eigenvectors)

        self.assertEqual(eigenvalues.Size(), 2)
        self.assertEqual(eigenvectors.Size1(), 2)
        self.assertEqual(eigenvectors.Size2(), 5)

        self.assertAlmostEqual(eigenvalues[0], 10.5, 7)
        self.assertAlmostEqual(eigenvalues[1], 12.0, 7)

        Kc = KratosMultiphysics.ComplexCompressedMatrix(K)
        Mc = KratosMultiphysics.ComplexCompressedMatrix(M)

        for i in range(eigenvalues.Size()):
            eigenvector = KratosMultiphysics.ComplexVector(n)
            for j in range(n):
                eigenvector[j] = eigenvectors[i,j]

            _aux_1 = KratosMultiphysics.ComplexVector(n,0)
            _aux_2 = KratosMultiphysics.ComplexVector(n,0)

            space.Mult(Kc,eigenvector,_aux_1)
            space.Mult(Mc,eigenvector,_aux_2)
            _aux_2 *= eigenvalues[i]

            self.assertVectorAlmostEqual(_aux_1, _aux_2)

    def test_real_general_gev_complex_result(self):

        space = KratosMultiphysics.UblasComplexSparseSpace()

        settings = KratosMultiphysics.Parameters('''{
            "solver_type" : "feast",
            "symmetric" : false,
            "subspace_size" : 2,
            "e_mid_re": 0.0,
            "e_mid_im": 0.0,
            "e_r": 30.0,
            "sort_eigenvalues": true,
            "sort_order": "li",
            "echo_level": 0
         }''')

        K = KratosMultiphysics.CompressedMatrix(2,2)
        K[0,0] = 1/sqrt(2)
        K[1,1] = 1

        n = K.Size1()

        M = KratosMultiphysics.CompressedMatrix(n, n)
        M[0,1] = 1
        M[1,0] = -1/sqrt(2)

        # create result containers
        eigenvalues = KratosMultiphysics.ComplexVector(n)
        eigenvectors = KratosMultiphysics.ComplexMatrix(n, 1)

        # Construct the solver
        eigen_solver = eigen_solver_factory.ConstructSolver(settings)

        # Solve
        eigen_solver.Solve(K, M, eigenvalues, eigenvectors)

        self.assertAlmostEqual(eigenvalues[0], 0.0+1.0j, 7)
        self.assertAlmostEqual(eigenvalues[1], 0.0-1.0j, 7)

        Kc = KratosMultiphysics.ComplexCompressedMatrix(K)
        Mc = KratosMultiphysics.ComplexCompressedMatrix(M)

        for i in range(eigenvalues.Size()):
            eigenvector = KratosMultiphysics.ComplexVector(n)
            for j in range(n):
                eigenvector[j] = eigenvectors[i,j]

            _aux_1 = KratosMultiphysics.ComplexVector(n,0)
            _aux_2 = KratosMultiphysics.ComplexVector(n,0)

            space.Mult(Kc,eigenvector,_aux_1)
            space.Mult(Mc,eigenvector,_aux_2)
            _aux_2 *= eigenvalues[i]

            self.assertVectorAlmostEqual(_aux_1, _aux_2)

    def test_complex_symmetric_gev(self):

        space = KratosMultiphysics.UblasComplexSparseSpace()

        settings = KratosMultiphysics.Parameters('''{
            "solver_type": "feast_complex",
            "symmetric": true,
            "number_of_eigenvalues": 3,
            "e_mid_re": 0.0,
            "e_mid_im": 0.0,
            "e_r": 0.16,
            "sort_eigenvalues": true,
            "sort_order": "lm",
            "echo_level": 0
         }''')

        K_real = KratosMultiphysics.CompressedMatrix()

        this_file_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(this_file_dir)))
        matrix_file_path = os.path.join(base_dir, "kratos", "tests", "test_files", "sparse_matrix_files", "A.mm")

        file_read = KratosMultiphysics.ReadMatrixMarketMatrix(matrix_file_path, K_real) # symmetric test matrix
        self.assertTrue(file_read, msg="The matrix file could not be read")

        K = KratosMultiphysics.ComplexCompressedMatrix(K_real)
        K *= .1j
        K += KratosMultiphysics.ComplexCompressedMatrix(K_real)

        n = K.Size1()
        self.assertEqual(n, 900)

        # create an identity matrix
        M = KratosMultiphysics.ComplexCompressedMatrix(n, n)
        for i in range(n):
            M[i, i] = 1.0

        # create result containers (they will be resized inside the solver)
        eigenvalues = KratosMultiphysics.ComplexVector(n)
        eigenvectors = KratosMultiphysics.ComplexMatrix(n, 1)

        # Construct the solver
        eigen_solver = eigen_solver_factory.ConstructSolver(settings)

        # Solve
        eigen_solver.Solve(K, M, eigenvalues, eigenvectors)

        for i in range(eigenvalues.Size()):
            eigenvector = KratosMultiphysics.ComplexVector(n)
            for j in range(n):
                eigenvector[j] = eigenvectors[i,j]

            _aux_1 = KratosMultiphysics.ComplexVector(n,0)
            _aux_2 = KratosMultiphysics.ComplexVector(n,0)

            space.Mult(K,eigenvector,_aux_1)
            space.Mult(M,eigenvector,_aux_2)
            _aux_2 *= eigenvalues[i]

            self.assertVectorAlmostEqual(_aux_1, _aux_2)

    def test_complex_general_gev(self):

        space = KratosMultiphysics.UblasComplexSparseSpace()

        settings = KratosMultiphysics.Parameters('''{
            "solver_type": "feast_complex",
            "symmetric": false,
            "subspace_size": 2,
            "e_mid_re": 10.0,
            "e_mid_im": 0.0,
            "e_r": 3.0,
            "sort_eigenvalues": true,
            "sort_order": "sm",
            "echo_level": 0
         }''')

        K_real = KratosMultiphysics.CompressedMatrix()

        this_file_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(this_file_dir)))
        matrix_file_path = os.path.join(base_dir, "kratos", "tests", "test_files", "sparse_matrix_files", "small_A.mm")

        file_read = KratosMultiphysics.ReadMatrixMarketMatrix(matrix_file_path, K_real) # symmetric test matrix
        self.assertTrue(file_read, msg="The matrix file could not be read")

        K = KratosMultiphysics.ComplexCompressedMatrix(K_real)
        K *= .1j
        K += KratosMultiphysics.ComplexCompressedMatrix(K_real)

        n = K.Size1()

        # create an identity matrix
        M = KratosMultiphysics.ComplexCompressedMatrix(n, n)
        for i in range(n):
            M[i, i] = 1.0

        # create result containers (they will be resized inside the solver)
        eigenvalues = KratosMultiphysics.ComplexVector(n)
        eigenvectors = KratosMultiphysics.ComplexMatrix(n, 1)

        # Construct the solver
        eigen_solver = eigen_solver_factory.ConstructSolver(settings)

        # Solve
        eigen_solver.Solve(K, M, eigenvalues, eigenvectors)

        self.assertEqual(eigenvalues.Size(), 2)
        self.assertEqual(eigenvectors.Size1(), 2)
        self.assertEqual(eigenvectors.Size2(), 5)

        self.assertAlmostEqual(eigenvalues[0], 10.5+1.05j, 7)
        self.assertAlmostEqual(eigenvalues[1], 12.0+1.2j, 7)

        for i in range(eigenvalues.Size()):
            eigenvector = KratosMultiphysics.ComplexVector(n)
            for j in range(n):
                eigenvector[j] = eigenvectors[i,j]

            _aux_1 = KratosMultiphysics.ComplexVector(n,0)
            _aux_2 = KratosMultiphysics.ComplexVector(n,0)

            space.Mult(K,eigenvector,_aux_1)
            space.Mult(M,eigenvector,_aux_2)
            _aux_2 *= eigenvalues[i]

            self.assertVectorAlmostEqual(_aux_1, _aux_2)


if __name__ == '__main__':
    KratosUnittest.main()