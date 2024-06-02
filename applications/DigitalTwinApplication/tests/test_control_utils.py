import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.DigitalTwinApplication as KratosDT

class TestControlUtils(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.n = 10
        cls.matrix = Kratos.Matrix(cls.n, cls.n, 0.0)
        cls.vector = Kratos.Vector(cls.n * (cls.n - 1) // 2 , 0.0)

        index = 0
        value = 1.0
        for i in range(cls.matrix.Size1()):
            for j in range(i + 1, cls.matrix.Size2()):
                cls.matrix[i, j] = value
                cls.matrix[j, i] = value
                cls.vector[index] = value
                value += 1.5
                index += 1

    def test_GetDistVectorSize(self):
        self.assertEqual(KratosDT.ControlUtils.GetDistVectorSize(1), 0)
        self.assertEqual(KratosDT.ControlUtils.GetDistVectorSize(2), 1)
        self.assertEqual(KratosDT.ControlUtils.GetDistVectorSize(3), 3)
        self.assertEqual(KratosDT.ControlUtils.GetDistVectorSize(4), 6)

    def test_GetDistIndexFromPairIndices(self):
        for i in range(self.n):
            for j in range(i + 1, self.n):
                dist_i = KratosDT.ControlUtils.GetDistIndexFromPairIndices(self.n, i, j)
                self.assertEqual(self.matrix[i, j], self.vector[dist_i])

    def test_GetPairIndicesFromDistIndex(self):
        m = KratosDT.ControlUtils.GetDistVectorSize(self.n)
        for dist_i in range(m):
            i, j = KratosDT.ControlUtils.GetPairIndicesFromDistIndex(self.n, dist_i)
            self.assertEqual(self.matrix[i, j], self.vector[dist_i])

if __name__ == '__main__':
    UnitTest.main()