# Import Kratos
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Additional imports
import os
import sys

import numpy as np

try:
    from scipy import sparse, io
    missing_scipy = False
except ImportError as e:
    missing_scipy = True

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestSparseMatrixSum(KratosUnittest.TestCase):

    def __sparse_matrix_sum(self, file_name = "test_files/sparse_matrix_files/A.mm"):
        # Read the matrices

        from scipy.io import mmread

        print('Reading matrix file ', file_name, file=sys.stderr)
        print('Reading matrix file ', file_name, file=sys.stdout)

        try:

            matrix_test = mmread(file_name)

        except Exception as e:

            print('Error reading matrix file ', file_name, file=sys.stderr)
            print('Error reading matrix file ', file_name, file=sys.stdout)

            raise e

        # A = KratosMultiphysics.CompressedMatrix()
        # B = KratosMultiphysics.CompressedMatrix()
        # KratosMultiphysics.ReadMatrixMarketMatrix(GetFilePath(file_name),A)
        # KratosMultiphysics.ReadMatrixMarketMatrix(GetFilePath(file_name),B)

        # Print A Matrix file
        with open(file_name, "r") as f:
            print(f.read(), file=sys.stderr)

        A_python = io.mmread(GetFilePath(file_name))
        A_python.toarray()
        B_python = io.mmread(GetFilePath(file_name))
        B_python.toarray()

        A_python = A_python + B_python

        # Solve
        KratosMultiphysics.SparseMatrixMultiplicationUtility.MatrixAdd(A, B, 1.0)

        for i, j in np.nditer(A_python.nonzero()):
            self.assertAlmostEqual(A[int(i), int(j)], A_python[int(i), int(j)])

    @KratosUnittest.skipIf(missing_scipy,"Missing python libraries (scipy)")
    def test_sparse_matrix_sum_small(self):
        self.__sparse_matrix_sum("test_files/sparse_matrix_files/small_A.mm")

    @KratosUnittest.skipIf(True,"Missing python libraries (scipy)")
    def test_sparse_matrix_sum_full(self):
        self.__sparse_matrix_sum()

class TestSparseMatrixTranspose(KratosUnittest.TestCase):

    def __sparse_matrix_transpose(self, file_name = "test_files/sparse_matrix_files/A.mm"):
        # Read the matrices
        A = KratosMultiphysics.CompressedMatrix()
        KratosMultiphysics.ReadMatrixMarketMatrix(GetFilePath(file_name),A)
        B = KratosMultiphysics.CompressedMatrix()

        A_python = io.mmread(GetFilePath(file_name))
        B_python = np.matrix.transpose(A_python.toarray())

        # Solve
        KratosMultiphysics.SparseMatrixMultiplicationUtility.TransposeMatrix(B, A, 1.0)

        for i, j in np.nditer(B_python.nonzero()):
            self.assertAlmostEqual(B[int(i), int(j)], A[int(j), int(i)])

    @KratosUnittest.skipIf(missing_scipy,"Missing python libraries (scipy)")
    def test_sparse_matrix_transpose_small(self):
        self.__sparse_matrix_transpose("test_files/sparse_matrix_files/small_A.mm")

    @KratosUnittest.skipIf(missing_scipy,"Missing python libraries (scipy)")
    def test_sparse_matrix_transpose_full(self):
        self.__sparse_matrix_transpose()

class TestSparseMatrixMultiplication(KratosUnittest.TestCase):

    def __sparse_matrix_multiplication(self, problem = "saad", file_name = "test_files/sparse_matrix_files/A.mm"):
        # Read the matrices
        A = KratosMultiphysics.CompressedMatrix()
        A2 = KratosMultiphysics.CompressedMatrix()
        KratosMultiphysics.ReadMatrixMarketMatrix(GetFilePath(file_name),A)

        A_python = io.mmread(GetFilePath(file_name))
        A_python.toarray()

        A2_python = np.dot(A_python, A_python)

        # Solve
        if problem == "saad":
            KratosMultiphysics.SparseMatrixMultiplicationUtility.MatrixMultiplicationSaad(A, A, A2)
        elif problem == "rmerge":
            KratosMultiphysics.SparseMatrixMultiplicationUtility.MatrixMultiplicationRMerge(A, A, A2)

        for i, j in np.nditer(A2_python.nonzero()):
            self.assertAlmostEqual(A2[int(i), int(j)], A2_python[int(i), int(j)], places=3)

    @KratosUnittest.skipIf(missing_scipy,"Missing python libraries (scipy)")
    def test_sparse_matrix_multiplication_saad_small(self):
        self.__sparse_matrix_multiplication("saad", "test_files/sparse_matrix_files/small_A.mm")

    @KratosUnittest.skipIf(missing_scipy,"Missing python libraries (scipy)")
    def test_sparse_matrix_multiplication_rmerge_small(self):
        self.__sparse_matrix_multiplication("rmerge", "test_files/sparse_matrix_files/small_A.mm")

    @KratosUnittest.skipIf(missing_scipy,"Missing python libraries (scipy)")
    def test_sparse_matrix_multiplication_saad_full(self):
        self.__sparse_matrix_multiplication("saad")

    @KratosUnittest.skipIf(missing_scipy,"Missing python libraries (scipy)")
    def test_sparse_matrix_multiplication_rmerge_full(self):
        self.__sparse_matrix_multiplication("rmerge")

if __name__ == '__main__':
    KratosUnittest.main()
