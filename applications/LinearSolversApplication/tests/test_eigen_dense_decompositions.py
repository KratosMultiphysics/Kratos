
from __future__ import print_function, absolute_import, division

import random

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.LinearSolversApplication as LinearSolversApplication

class TestEigenDenseDecompositions(KratosUnittest.TestCase):
    def _square_matrix_3x3(self):
        A = KratosMultiphysics.Matrix(3,3)
        A[0,0] = 0.57690
        A[0,1] = 0.28760
        A[0,2] = 0.63942
        A[1,0] = 0.72886
        A[1,1] = 0.40541
        A[1,2] = 0.13415
        A[2,0] = 0.81972
        A[2,1] = 0.54501
        A[2,2] = 0.28974

        return A

    def _random_matrix_mxn(self, m, n):
        A = KratosMultiphysics.Matrix(m,n)
        for i in range(m):
            for j in range(n):
                A[i,j] = random.random()

        return A

    def _execute_eigen_dense_decomposition_test(self, decomposition, generate_matrix_function, thin_u_v = False):
        # Generate the sampling matrix
        A = generate_matrix_function()

        # Compute the SVD decomposition
        S = KratosMultiphysics.Vector()
        U = KratosMultiphysics.Matrix()
        V = KratosMultiphysics.Matrix()
        settings = KratosMultiphysics.Parameters('''{}''')
        if thin_u_v:
            settings.AddBool("compute_thin_u", True)
            settings.AddBool("compute_thin_v", True)
        decomposition.Compute(A, S, U, V, settings)

        # Reconstruct the sampling matrix with the SVD arrays
        aux_S_mat = KratosMultiphysics.Matrix(S.Size(), V.Size2(), 0.0)
        for i in range(S.Size()):
            aux_S_mat[i,i] = S[i]

        m = A.Size1()
        n = A.Size2()
        A_svd = KratosMultiphysics.Matrix(m, n, 0.0)
        for i in range(m):
            for j in range(n):
                for k in range(U.Size2()):
                    for m in range(aux_S_mat.Size2()):
                        A_svd[i,j] += U[i,k] * aux_S_mat[k,m] * V[j,m]

        # Check that input matrix equals the decomposition reconstructed one
        for i in range(A.Size1()):
            for j in range(A.Size2()):
                self.assertAlmostEqual(A[i,j], A_svd[i,j], 7)

    def test_eigen_dense_BDC_SVD_3x3(self):
        decomposition = LinearSolversApplication.EigenDenseBDCSVD()
        self._execute_eigen_dense_decomposition_test(decomposition, self._square_matrix_3x3)

    def test_eigen_dense_BDC_SVD_5x4_random(self):
        decomposition = LinearSolversApplication.EigenDenseBDCSVD()
        self._execute_eigen_dense_decomposition_test(decomposition, lambda : self._random_matrix_mxn(5,4))

    def test_eigen_dense_BDC_SVD_3x10_random(self):
        decomposition = LinearSolversApplication.EigenDenseBDCSVD()
        self._execute_eigen_dense_decomposition_test(decomposition, lambda : self._random_matrix_mxn(3,10))

    def test_eigen_dense_BDC_SVD_5x4_random_thin(self):
        decomposition = LinearSolversApplication.EigenDenseBDCSVD()
        self._execute_eigen_dense_decomposition_test(decomposition, lambda : self._random_matrix_mxn(5,4), True)

    def test_eigen_dense_BDC_SVD_3x10_random_thin(self):
        decomposition = LinearSolversApplication.EigenDenseBDCSVD()
        self._execute_eigen_dense_decomposition_test(decomposition, lambda : self._random_matrix_mxn(3,10), True)

    def test_eigen_dense_Jacobi_SVD_3x3(self):
        decomposition = LinearSolversApplication.EigenDenseJacobiSVD()
        self._execute_eigen_dense_decomposition_test(decomposition, self._square_matrix_3x3)

    def test_eigen_dense_Jacobi_SVD_5x4_random(self):
        decomposition = LinearSolversApplication.EigenDenseJacobiSVD()
        self._execute_eigen_dense_decomposition_test(decomposition, lambda : self._random_matrix_mxn(5,4))

    def test_eigen_dense_Jacobi_SVD_3x10_random(self):
        decomposition = LinearSolversApplication.EigenDenseJacobiSVD()
        self._execute_eigen_dense_decomposition_test(decomposition, lambda : self._random_matrix_mxn(3,10))

    def test_eigen_dense_Jacobi_SVD_5x4_random_thin(self):
        decomposition = LinearSolversApplication.EigenDenseJacobiSVD()
        self._execute_eigen_dense_decomposition_test(decomposition, lambda : self._random_matrix_mxn(5,4), True)

    def test_eigen_dense_Jacobi_SVD_3x10_random_thin(self):
        decomposition = LinearSolversApplication.EigenDenseJacobiSVD()
        self._execute_eigen_dense_decomposition_test(decomposition, lambda : self._random_matrix_mxn(3,10), True)

if __name__ == '__main__':
    KratosUnittest.main()
