import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

class TestMatrixMarketInterface(KratosUnittest.TestCase):

    def test_mm_matrix_io(self):
        A = KratosMultiphysics.CompressedMatrix(3,3)
        A[0,0] = 1.0
        A[1,1] = 2.0
        A[2,2] = 3.0
        A[0,2] = 4.0

        KratosMultiphysics.WriteMatrixMarketMatrix('A.mm', A, False)

        B = KratosMultiphysics.CompressedMatrix()
        KratosMultiphysics.ReadMatrixMarketMatrix('A.mm', B)

        self.assertMatrixAlmostEqual(A, B)
        kratos_utils.DeleteFileIfExisting('A.mm')

    def test_mm_vector_io(self):
        a = KratosMultiphysics.Vector(5,0)
        a[0] = 1.0
        a[3] = 5.0

        KratosMultiphysics.WriteMatrixMarketVector('a.mm', a)

        b = KratosMultiphysics.Vector()
        KratosMultiphysics.ReadMatrixMarketVector('a.mm', b)

        self.assertVectorAlmostEqual(a, b)
        kratos_utils.DeleteFileIfExisting("a.mm")

    def test_mm_matrix_io_cplx(self):
        A = KratosMultiphysics.ComplexCompressedMatrix(3,3)
        A[0,0] = 1.0-1.0j
        A[1,1] = 2.0+2.0j
        A[2,2] = 3.0-3.0j
        A[0,2] = 4.0+4.0j

        KratosMultiphysics.WriteMatrixMarketMatrix('A.mm', A, False)

        B = KratosMultiphysics.ComplexCompressedMatrix()
        KratosMultiphysics.ReadMatrixMarketMatrix('A.mm', B)

        self.assertMatrixAlmostEqual(A, B)
        kratos_utils.DeleteFileIfExisting('A.mm')

    def test_mm_vector_io_cplx(self):
        a = KratosMultiphysics.ComplexVector(5,0)
        a[0] = 1.0-1.0j
        a[3] = 5.0+3.0j

        KratosMultiphysics.WriteMatrixMarketVector('a.mm', a)

        b = KratosMultiphysics.ComplexVector()
        KratosMultiphysics.ReadMatrixMarketVector('a.mm', b)

        self.assertVectorAlmostEqual(a, b)
        kratos_utils.DeleteFileIfExisting("a.mm")

if __name__ == '__main__':
    KratosUnittest.main()
