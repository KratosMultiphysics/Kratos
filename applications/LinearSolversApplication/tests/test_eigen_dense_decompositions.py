
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

    def _execute_eigen_dense_svd_test(self, decomposition, generate_matrix_function, thin_u_v = False):
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

    def _execute_eigen_dense_qr_decomposition_test(self, decomposition, generate_matrix_function, check_q_times_r):
        # Generate the sampling matrix
        A = generate_matrix_function()
        A_rows = A.Size1()
        A_cols = A.Size2()

        # Compute the decomposition
        decomposition.Compute(A)

        # Solve with the input matrix as RHS
        aux_sol = KratosMultiphysics.Matrix()
        decomposition.Solve(A, aux_sol)

        # Check that the obtained solution equals the identity
        m = aux_sol.Size1()
        n = aux_sol.Size2()
        for i in range(m):
            for j in range(n):
                ref_val = 1.0 if ((i == j) and (i < A_cols)) else 0.0
                self.assertAlmostEqual(aux_sol[i,j], ref_val, 7)

        # Check the QR matrix reconstruction
        # Note that we check AP=QR since R is only available for the pivoting QR decompositions
        if check_q_times_r:
            Q = KratosMultiphysics.Matrix()
            R = KratosMultiphysics.Matrix()
            P = KratosMultiphysics.Matrix()
            decomposition.MatrixQ(Q)
            decomposition.MatrixR(R)
            decomposition.MatrixP(P)
            A_decomp = KratosMultiphysics.Matrix(A_rows,A_cols,0.0)
            A_permut = KratosMultiphysics.Matrix(A_rows,A_cols,0.0)

            for i in range(A_rows):
                for j in range(A_cols):
                    for k in range(R.Size1()):
                        A_decomp[i,j] += Q[i,k]*R[k,j]

            for i in range(A_rows):
                for j in range(A_cols):
                    for k in range(P.Size2()):
                        A_permut[i,j] += A[i,k]*P[k,j]

            for i in range(m):
                for j in range(n):
                    self.assertAlmostEqual(A_decomp[i,j], A_permut[i,j], 7)

    def test_eigen_dense_BDC_SVD_3x3(self):
        decomposition = LinearSolversApplication.EigenDenseBDCSVD()
        self._execute_eigen_dense_svd_test(decomposition, self._square_matrix_3x3)

    def test_eigen_dense_BDC_SVD_5x4_random(self):
        decomposition = LinearSolversApplication.EigenDenseBDCSVD()
        self._execute_eigen_dense_svd_test(decomposition, lambda : self._random_matrix_mxn(5,4))

    def test_eigen_dense_BDC_SVD_3x10_random(self):
        decomposition = LinearSolversApplication.EigenDenseBDCSVD()
        self._execute_eigen_dense_svd_test(decomposition, lambda : self._random_matrix_mxn(3,10))

    def test_eigen_dense_BDC_SVD_5x4_random_thin(self):
        decomposition = LinearSolversApplication.EigenDenseBDCSVD()
        self._execute_eigen_dense_svd_test(decomposition, lambda : self._random_matrix_mxn(5,4), True)

    def test_eigen_dense_BDC_SVD_3x10_random_thin(self):
        decomposition = LinearSolversApplication.EigenDenseBDCSVD()
        self._execute_eigen_dense_svd_test(decomposition, lambda : self._random_matrix_mxn(3,10), True)

    def test_eigen_dense_Jacobi_SVD_3x3(self):
        decomposition = LinearSolversApplication.EigenDenseJacobiSVD()
        self._execute_eigen_dense_svd_test(decomposition, self._square_matrix_3x3)

    def test_eigen_dense_Jacobi_SVD_5x4_random(self):
        decomposition = LinearSolversApplication.EigenDenseJacobiSVD()
        self._execute_eigen_dense_svd_test(decomposition, lambda : self._random_matrix_mxn(5,4))

    def test_eigen_dense_Jacobi_SVD_3x10_random(self):
        decomposition = LinearSolversApplication.EigenDenseJacobiSVD()
        self._execute_eigen_dense_svd_test(decomposition, lambda : self._random_matrix_mxn(3,10))

    def test_eigen_dense_Jacobi_SVD_5x4_random_thin(self):
        decomposition = LinearSolversApplication.EigenDenseJacobiSVD()
        self._execute_eigen_dense_svd_test(decomposition, lambda : self._random_matrix_mxn(5,4), True)

    def test_eigen_dense_Jacobi_SVD_3x10_random_thin(self):
        decomposition = LinearSolversApplication.EigenDenseJacobiSVD()
        self._execute_eigen_dense_svd_test(decomposition, lambda : self._random_matrix_mxn(3,10), True)

    def test_eigen_dense_Householder_QR_decomposition_3x3(self):
        decomposition = LinearSolversApplication.EigenDenseHouseholderQRDecomposition()
        self._execute_eigen_dense_qr_decomposition_test(decomposition, self._square_matrix_3x3, False)

    def test_eigen_dense_Householder_QR_decomposition_5x4_random(self):
        decomposition = LinearSolversApplication.EigenDenseHouseholderQRDecomposition()
        self._execute_eigen_dense_qr_decomposition_test(decomposition, lambda : self._random_matrix_mxn(5,4), False)

    def test_eigen_dense_Householder_QR_decomposition_10x3_random(self):
        decomposition = LinearSolversApplication.EigenDenseHouseholderQRDecomposition()
        self._execute_eigen_dense_qr_decomposition_test(decomposition, lambda : self._random_matrix_mxn(10,3), False)

    def test_eigen_dense_column_pivoting_Householder_QR_decomposition_3x3(self):
        decomposition = LinearSolversApplication.EigenDenseColumnPivotingHouseholderQRDecomposition()
        self._execute_eigen_dense_qr_decomposition_test(decomposition, self._square_matrix_3x3, True)

    def test_eigen_dense_column_pivoting_Householder_QR_decomposition_5x4_random(self):
        decomposition = LinearSolversApplication.EigenDenseColumnPivotingHouseholderQRDecomposition()
        self._execute_eigen_dense_qr_decomposition_test(decomposition, lambda : self._random_matrix_mxn(5,4), True)

    def test_eigen_dense_column_pivoting_Householder_QR_decomposition_10x3_random(self):
        decomposition = LinearSolversApplication.EigenDenseColumnPivotingHouseholderQRDecomposition()
        self._execute_eigen_dense_qr_decomposition_test(decomposition, lambda : self._random_matrix_mxn(10,3), True)

if __name__ == '__main__':
    KratosUnittest.main()
