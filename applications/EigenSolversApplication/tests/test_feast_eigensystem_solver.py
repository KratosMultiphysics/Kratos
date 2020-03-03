from __future__ import print_function, absolute_import, division
import os
import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import eigen_solver_factory
import KratosMultiphysics.EigenSolversApplication as EigenSolversApplication   

class TestFeastEigensystemSolver(KratosUnittest.TestCase):
    @KratosUnittest.skipUnless(hasattr(EigenSolversApplication,'FEASTEigensystemSolver'),"FEAST not found, skipping.")
    def test_real_symmetric_gev(self):

        space = KratosMultiphysics.UblasSparseSpace()

        settings = KratosMultiphysics.Parameters('''{
            "solver_type": "eigen_feast",
            "number_of_eigenvalues": 3,
            "search_lowest_eigenvalues": true,
            "e_min" : 0.0,
            "e_max" : 0.2,
            "echo_level": 0
         }''')

        K = KratosMultiphysics.CompressedMatrix()

        this_file_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(this_file_dir)))
        matrix_file_path = os.path.join(base_dir, "kratos", "tests", "auxiliar_files_for_python_unittest", "sparse_matrix_files", "A.mm")

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
                eigenvector[j] = eigenvectors[j,i]

            _aux_1 = KratosMultiphysics.Vector(n,0)
            _aux_2 = KratosMultiphysics.Vector(n,0)

            space.Mult(K,eigenvector,_aux_1)
            space.Mult(M,eigenvector,_aux_2)
            _aux_2 *= eigenvalues[i]

            self.assertVectorAlmostEqual(_aux_1, _aux_2)

    @KratosUnittest.skipUnless(hasattr(EigenSolversApplication,'ComplexFEASTEigensystemSolver'),"FEAST not found, skipping.")
    def test_complex_symmetric_gev(self):

        space = KratosMultiphysics.UblasComplexSparseSpace()

        settings = KratosMultiphysics.Parameters('''{
            "solver_type": "eigen_feast_complex",
            "number_of_eigenvalues": 3,
            "e_mid_re" : 0.0,
            "e_mid_im" : 0.0,
            "e_r" : 0.16,
            "echo_level": 0
         }''')

        K_real = KratosMultiphysics.CompressedMatrix()

        this_file_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(this_file_dir)))
        matrix_file_path = os.path.join(base_dir, "kratos", "tests", "auxiliar_files_for_python_unittest", "sparse_matrix_files", "A.mm")

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
                eigenvector[j] = eigenvectors[j,i]

            _aux_1 = KratosMultiphysics.ComplexVector(n,0)
            _aux_2 = KratosMultiphysics.ComplexVector(n,0)

            space.Mult(K,eigenvector,_aux_1)
            space.Mult(M,eigenvector,_aux_2)
            _aux_2 *= eigenvalues[i]

            self.assertVectorAlmostEqual(_aux_1, _aux_2)


if __name__ == '__main__':
    KratosUnittest.main()