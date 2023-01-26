import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM

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

    def test_list_of_list_construction_error_hand(self):
        with self.assertRaisesRegex(RuntimeError, r'Error: Wrong size of a row 1! Expected 2, got 3'):
            KM.Matrix([[1, 2], [4, 5, 6]])


if __name__ == '__main__':
    KratosUnittest.main()
