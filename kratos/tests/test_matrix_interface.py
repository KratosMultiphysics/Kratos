from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *
import math


class TestMatrixInterface(KratosUnittest.TestCase):
    def test_assignement(self):
        a = Matrix(2,3)

        self.assertEqual(2,a.Size1())
        self.assertEqual(3,a.Size2())

        for i in range(a.Size1()):
            for j in range(a.Size2()):
                a[i,j] = i+j

        for j in range(a.Size2()):
            for i in range(a.Size1()):
                self.assertEqual(a[i,j], i+j)

    def test_matrix_vector(self):
        A = Matrix(4,3)

        for i in range(A.Size1()):
            for j in range(A.Size2()):
                A[i,j] = i

        #matrix vector
        b = Vector(3,1.0)
        c = A*b
        for i in range(len(c)):
            self.assertEqual(c[i],i*A.Size2())

        #matrix array_1d<double,3>
        b = Array3(1.0)
        c = A*b
        for i in range(len(c)):
            self.assertEqual(c[i],i*A.Size2())


    def test_matrix_sum(self):
        A = Matrix(2,3,1.0)
        B = Matrix(2,3,2.0)
        C = A+B

        for i in range(A.Size1()):
            for j in range(A.Size2()):
                self.assertEqual(C[i,j], 3.0)

        A += B

        for i in range(A.Size1()):
            for j in range(A.Size2()):
                self.assertEqual(C[i,j], A[i,j])

    def test_matrix_diff(self):
        A = Matrix(2,3,1.0)
        B = Matrix(2,3,2.0)
        C = A-B

        for i in range(A.Size1()):
            for j in range(A.Size2()):
                self.assertEqual(C[i,j], -1.0)

        A -= B
        for i in range(A.Size1()):
            for j in range(A.Size2()):
                self.assertEqual(C[i,j], A[i,j])

    def test_scalar_prod(self):
        A = Matrix(2,3,2.0)
        C = A*2.0

        for i in range(A.Size1()):
            for j in range(A.Size2()):
                self.assertEqual(C[i,j], 4.0)

        A *= 2.0
        print(A)
        for i in range(A.Size1()):
            for j in range(A.Size2()):
                self.assertEqual(C[i,j], A[i,j])



if __name__ == '__main__':
    KratosUnittest.main()
