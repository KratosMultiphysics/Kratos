import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestDynamicDimensionalArray(KratosUnittest.TestCase):
    def test_DynamicDimensionalArrayVanishingDims(self):
        kratos_array = Kratos.DoubleDynamicDimensionalArray([10, 0, 5])
        numpy_array = kratos_array.to_numpy()
        self.assertEqual(numpy_array.shape, (10, 0, 5))
        self.assertEqual(numpy_array.size, 0)

    def test_DynamicDimensionalArray0(self):
        kratos_array = Kratos.DoubleDynamicDimensionalArray([0, 0])
        numpy_array = kratos_array.to_numpy()
        self.assertEqual(numpy_array.shape, (0, 0))
        self.assertEqual(numpy_array.size, 0)

    def test_DynamicDimensionalArray1(self):
        np_array = np.ones((10,3))
        kratos_array = Kratos.DoubleDynamicDimensionalArray(np_array)
        del kratos_array # should not delete the numpy array data

        self.assertEqual(np.sum(np_array), 30)

    def test_DynamicDimensionalArray2(self):
        kratos_array = Kratos.DoubleDynamicDimensionalArray([10, 2, 2])

        # assign value
        numpy_array = kratos_array.to_numpy()
        numpy_array[:] = np.ones((10, 2, 2))

        numpy_array_2 = kratos_array.to_numpy()
        self.assertEqual(np.linalg.norm(numpy_array_2 - numpy_array), 0.0)

        # following will also keep the capsule of numpy_array_2 within the lifetime of numpy_array_3
        # hence the internal data of the kratos_array will be kept alive.
        numpy_array_3 = np.array(numpy_array_2, copy=False)
        numpy_array_4 = kratos_array.to_numpy()

        del kratos_array # delete the kratos array. But still the numpy array
                         # should be there, since the internal data is shared
                         # with memory management.

        self.assertEqual(np.linalg.norm(numpy_array_2 - numpy_array), 0.0)
        self.assertEqual(np.linalg.norm(numpy_array_3 - numpy_array_4), 0.0)

        numpy_array[0,0,0] = -2.0
        self.assertEqual(np.sum(numpy_array), 37)
        self.assertEqual(np.sum(numpy_array_2), 37)
        self.assertEqual(np.sum(numpy_array_3), 37)
        self.assertEqual(np.sum(numpy_array_4), 37)

    def test_DynamicDimensionalArray3(self):
        np_array = np.ones((10,3))
        kratos_array = Kratos.DoubleDynamicDimensionalArray(np_array, copy=False)

        np_array[5,1] = -20
        self.assertEqual(np.sum(kratos_array.to_numpy()), 9)

    def test_DynamicDimensionalArray4(self):
        for dtype in [bool, np.int8, np.int16, np.int32, np.int64,  np.float16, np.float32, np.uint8, np.uint16, np.uint32, np.uint64, np.uint]:
            with self.assertRaises(RuntimeError):
                np_array = np.ones((2,2), dtype=dtype)
                Kratos.DoubleDynamicDimensionalArray(np_array, copy=False)

        np_array = np.ones((5,5), dtype=np.float64)
        k_array = Kratos.DoubleDynamicDimensionalArray(np_array, copy=False)
        np_array[2,3] = 10
        self.assertEqual(k_array.to_numpy()[2,3], 10)

    def test_DynamicDimensionalArray5(self):
        for dtype in [bool, np.int8, np.int16, np.int32, np.int64,  np.float16, np.float32, np.float64, np.uint16, np.uint32, np.uint64, np.uint]:
            with self.assertRaises(RuntimeError):
                np_array = np.ones((2,2), dtype=dtype)
                Kratos.UIntDynamicDimensionalArray(np_array, copy=False)

        np_array = np.ones((5,5), dtype=np.uint8)
        k_array = Kratos.UIntDynamicDimensionalArray(np_array, copy=False)
        np_array[2,3] = 10
        self.assertEqual(k_array.to_numpy()[2,3], 10)

    def test_DynamicDimensionalArray6(self):
        for dtype in [bool, np.int8, np.int16, np.int64,  np.float16, np.float32, np.float64, np.uint8, np.uint16, np.uint32, np.uint64, np.uint]:
            with self.assertRaises(RuntimeError):
                np_array = np.ones((5,5), dtype=dtype)
                Kratos.IntDynamicDimensionalArray(np_array, copy=False)

        np_array = np.ones((5,5), dtype=np.int32)
        k_array = Kratos.IntDynamicDimensionalArray(np_array, copy=False)
        np_array[2,3] = 10
        self.assertEqual(k_array.to_numpy()[2,3], 10)

    def test_DynamicDimensionalArray6(self):
        for dtype in [np.int8, np.int16, np.int32, np.int64,  np.float16, np.float32, np.float64, np.uint8, np.uint16, np.uint32, np.uint64, np.uint]:
            with self.assertRaises(RuntimeError):
                np_array = np.ones((5,5), dtype=dtype)
                Kratos.BoolDynamicDimensionalArray(np_array, copy=False)

        np_array = np.ones((5,5), dtype=bool)
        k_array = Kratos.BoolDynamicDimensionalArray(np_array, copy=False)
        np_array[2,3] = False
        self.assertEqual(k_array.to_numpy()[2,3], False)

if __name__ == '__main__':
    KratosUnittest.main()