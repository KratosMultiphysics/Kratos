import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM

import numpy

class TestMatrixInterface(KratosUnittest.TestCase):

    def test_len(self):
        # creation
        A = KM.Matrix(9, 8)
        self.assertEqual(72, len(A))

        # after resizing
        A.Resize(5, 1)
        self.assertEqual(5, len(A))

    def test_copy(self):
        A = KM.Matrix(2,3)
        A.fill(1.0)

        b = KM.Vector(3)
        b[0] = 1.0
        b[1] = 2.0
        b[2] = 3.0

        A2 = KM.Matrix(A)
        b2 = KM.Vector(b)

        self.assertVectorAlmostEqual(b, b2)
        self.assertMatrixAlmostEqual(A, A2)

        A[0,1] = 2.0
        self.assertEqual(1.0, A2[0,1])
        b[0] = 2.0
        self.assertEqual(1.0, b2[0])

    def test_assignement(self):
        A = KM.Matrix(2,3)

        self.assertEqual(2,A.Size1())
        self.assertEqual(3,A.Size2())

        for i in range(A.Size1()):
            for j in range(A.Size2()):
                A[i,j] = i+j

        for j in range(A.Size2()):
            for i in range(A.Size1()):
                self.assertEqual(A[i,j], i+j)

    def test_matrix_vector(self):
        A = KM.Matrix(4,3)

        for i in range(A.Size1()):
            for j in range(A.Size2()):
                A[i,j] = i

        #matrix vector
        b = KM.Vector(3)
        b.fill(1.0)
        c = A*b
        for i in range(len(c)):
            self.assertEqual(c[i],i*A.Size2())

        #matrix array_1d<double,3>
        b = KM.Array3(1.0)
        c = A*b
        for i in range(len(c)):
            self.assertEqual(c[i],i*A.Size2())


    def test_matrix_sum(self):
        A = KM.Matrix(2,3)
        A.fill(1.0)
        B = KM.Matrix(2,3)
        B.fill(2.0)
        C = A+B

        for i in range(A.Size1()):
            for j in range(A.Size2()):
                self.assertEqual(C[i,j], 3.0)

        A += B
        self.assertMatrixAlmostEqual(C, A)

    def test_matrix_diff(self):
        A = KM.Matrix(2,3)
        A.fill(1.0)
        B = KM.Matrix(2,3)
        B.fill(2.0)
        C = A-B

        for i in range(A.Size1()):
            for j in range(A.Size2()):
                self.assertEqual(C[i,j], -1.0)

        A -= B
        self.assertMatrixAlmostEqual(C, A)

    def test_scalar_prod(self):
        A = KM.Matrix(2,3)
        A.fill(2.0)
        C = A*2.0

        for i in range(A.Size1()):
            for j in range(A.Size2()):
                self.assertEqual(C[i,j], 4.0)

        A *= 2.0
        self.assertMatrixAlmostEqual(C, A)

    def test_matrix_transpose(self):
        A = KM.Matrix(2,3)
        B = KM.Matrix(3,2)
        for i in range(A.Size1()):
            for j in range(A.Size2()):
                A[i,j] = i
        B = A.transpose()

        for i in range(A.Size1()):
            for j in range(A.Size2()):
                self.assertEqual(A[i,j], B[j,i])

    def test_empty_list_construction(self):
        A = KM.Matrix([])
        self.assertMatrixAlmostEqual(A, KM.Matrix(0,0))

    def test_empty_list_of_list_construction(self):
        B = KM.Matrix([[]])
        self.assertMatrixAlmostEqual(B, KM.Matrix(0,0))

    def test_rectangular_list_of_list_construction(self):
        C = KM.Matrix([[2, 2, 2, 2]])
        self.assertEqual(4, len(C))
        self.assertEqual(1, C.Size1())
        self.assertEqual(4, C.Size2())

        for i in range(C.Size1()):
            for j in range(C.Size2()):
                self.assertEqual(C[i,j], 2)

    def test_sqaure_list_of_list_construction(self):
        l2 = [[0, 0.1, 0.2], [1, 1.1, 1.2], [2, 2.1, 2.2]]
        D = KM.Matrix(l2)
        self.assertEqual(9, len(D))
        self.assertEqual(3, D.Size1())
        self.assertEqual(3, D.Size2())

        for i in range(D.Size1()):
            for j in range(D.Size2()):
                self.assertEqual(D[i,j], i + (j*0.1))

    def test_buffer_protocol(self) -> None:
        # + ------------ +
        # |  0  1  2  3  |
        # |  4  5  6  7  |
        # |  8  9 10 11  |
        # | 12 13 14 15  |
        # | 16 17 18 19  |
        # + ------------ +
        numpy_matrix = numpy.arange(20.0, dtype=numpy.float64).reshape((5, 4))

        # Test default casting
        kratos_matrix = KM.Matrix(numpy_matrix)
        reference_matrix = KM.Matrix([[ 0.0,  1.0,  2.0,  3.0],
                                      [ 4.0,  5.0,  6.0,  7.0],
                                      [ 8.0,  9.0, 10.0, 11.0],
                                      [12.0, 13.0, 14.0, 15.0],
                                      [16.0, 17.0, 18.0, 19.0]])
        self.assertMatrixAlmostEqual(kratos_matrix, reference_matrix)

        # Test slicing
        # + ------ +
        # |  5  6  |
        # |  9 10  |
        # | 13 14  |
        # + ------ +
        sliced_numpy_matrix = numpy_matrix[1:4,1:3]
        kratos_matrix = KM.Matrix(sliced_numpy_matrix)
        reference_matrix = KM.Matrix([[ 5.0,  6.0],
                                      [ 9.0, 10.0],
                                      [13.0, 14.0]])
        self.assertMatrixAlmostEqual(kratos_matrix, reference_matrix)

        # Test slicing with steps
        # + ----- +
        # |  1  3 |
        # | 13 15 |
        # + ----- +
        sliced_numpy_matrix = numpy_matrix[0::3, 1::2]
        kratos_matrix = KM.Matrix(sliced_numpy_matrix)
        reference_matrix = KM.Matrix([[ 1.0,  3.0],
                                      [13.0, 15.0]])
        self.assertMatrixAlmostEqual(kratos_matrix, reference_matrix)

        # Test empty slice
        # ++
        # ++
        sliced_numpy_matrix = numpy_matrix[10:10,20:20]
        kratos_matrix = KM.Matrix(sliced_numpy_matrix)
        reference_matrix = KM.Matrix([])
        self.assertMatrixAlmostEqual(kratos_matrix, reference_matrix)

    def test_list_of_list_construction_error_hand(self):
        with self.assertRaisesRegex(RuntimeError, r'Error: Wrong size of a row 1! Expected 2, got 3'):
            KM.Matrix([[1, 2], [4, 5, 6]])


if __name__ == '__main__':
    KratosUnittest.main()
