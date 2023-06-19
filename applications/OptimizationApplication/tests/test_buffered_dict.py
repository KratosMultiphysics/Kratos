
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict

class TestBufferedDict(kratos_unittest.TestCase):
    def test_AdvanceStep(self):
        buffered_data = BufferedDict(3)

        buffered_data.SetValue("step", 1)
        self.assertEqual(buffered_data.GetValue("step"), 1)

        buffered_data.AdvanceStep()
        buffered_data.SetValue("step", 2)
        self.assertEqual(buffered_data.GetValue("step"), 2)
        self.assertEqual(buffered_data.GetValue("step", 1), 1)

        buffered_data.AdvanceStep()
        buffered_data.SetValue("step", 3)
        self.assertEqual(buffered_data.GetValue("step"), 3)
        self.assertEqual(buffered_data.GetValue("step", 1), 2)
        self.assertEqual(buffered_data.GetValue("step", 2), 1)

        buffered_data.AdvanceStep()
        buffered_data.SetValue("step", 4)
        self.assertEqual(buffered_data.GetValue("step"), 4)
        self.assertEqual(buffered_data.GetValue("step", 1), 3)
        self.assertEqual(buffered_data.GetValue("step", 2), 2)
        with self.assertRaises(RuntimeError):
            _ = buffered_data.GetValue("step", 3)

    def test_AdvanceStepSet(self):
        buffered_data = BufferedDict(3)
        buffered_data.SetValue("test_1/test_sub_1", BufferedDict(2))

        self.assertEqual(buffered_data.GetValue("test_1").GetBufferSize(), 3)
        self.assertEqual(buffered_data.GetValue("test_1/test_sub_1").GetBufferSize(), 2)

        buffered_data.SetValue("test_1/int", 1)
        buffered_data.SetValue("test_1/test_sub_1/int", 2)

        buffered_data.AdvanceStep()
        buffered_data.SetValue("test_1/int", 3)
        buffered_data.SetValue("test_1/test_sub_1/int", 4)

        buffered_data.AdvanceStep()
        buffered_data.SetValue("test_1/int", 5)
        buffered_data.SetValue("test_1/test_sub_1/int", 6)

        buffered_data.AdvanceStep()
        buffered_data.SetValue("test_1/int", 7)
        buffered_data.SetValue("test_1/test_sub_1/int", 8)

        buffered_data.AdvanceStep()
        buffered_data.SetValue("test_1/int", 9)
        buffered_data.SetValue("test_1/test_sub_1/int", 10)

        buffered_data.AdvanceStep()
        buffered_data.SetValue("test_1/int", 11)
        buffered_data.SetValue("test_1/test_sub_1/int", 12)

        self.assertEqual(buffered_data.GetValue("test_1/int", 0), 11)
        self.assertEqual(buffered_data.GetValue("test_1/int", 1), 9)
        self.assertEqual(buffered_data.GetValue("test_1/int", 2), 7)
        with self.assertRaises(RuntimeError):
            _ = buffered_data.GetValue("test_1/int", 3)

        self.assertEqual(buffered_data.GetValue("test_1/test_sub_1/int", 0), 12)
        self.assertEqual(buffered_data.GetValue("test_1/test_sub_1/int", 1), 10)
        with self.assertRaises(RuntimeError):
            _ = buffered_data.GetValue("test_1/test_sub_1/int", 3)

    def test_HasValue(self):
        buffered_data = BufferedDict(2)
        buffered_data.SetValue("data/sub_1", 1)
        buffered_data.SetValue("data/sub_2", 2, 1)
        buffered_data.SetValue("data/sub_3/sub_sub1", 3, 1)

        self.assertTrue(buffered_data.HasValue("data/sub_2", 1))
        self.assertTrue(buffered_data.GetValue("data").HasValue("sub_3/sub_sub1", 1))
        self.assertTrue(buffered_data.HasValue("data/sub_2", 1))
        self.assertTrue(buffered_data.HasValue("data/sub_3", 0))
        self.assertTrue(buffered_data.HasValue("data/sub_3", 1))
        self.assertFalse(buffered_data.HasValue("data/sub_2"))
        self.assertFalse(buffered_data.HasValue("data_2/sub_2"))

    def test_RemoveValue(self):
        buffered_data = BufferedDict(3)
        buffered_data.SetValue("data/sub_1", 1)
        buffered_data.SetValue("data/sub_1", 2, 1)
        buffered_data.SetValue("data/sub_1", 3, 2)

        self.assertTrue(buffered_data.HasValue("data/sub_1", 0))
        self.assertTrue(buffered_data.HasValue("data/sub_1", 1))
        self.assertTrue(buffered_data.HasValue("data/sub_1", 2))

        buffered_data.RemoveValue("data/sub_1", 0)
        self.assertFalse(buffered_data.HasValue("data/sub_1", 0))
        self.assertTrue(buffered_data.HasValue("data/sub_1", 1))
        self.assertTrue(buffered_data.HasValue("data/sub_1", 2))

        buffered_data.RemoveValue("data/sub_1", 1)
        self.assertFalse(buffered_data.HasValue("data/sub_1", 0))
        self.assertFalse(buffered_data.HasValue("data/sub_1", 1))
        self.assertTrue(buffered_data.HasValue("data/sub_1", 2))

        buffered_data.RemoveValue("data/sub_1", 2)
        self.assertFalse(buffered_data.HasValue("data/sub_1", 0))
        self.assertFalse(buffered_data.HasValue("data/sub_1", 1))
        self.assertFalse(buffered_data.HasValue("data/sub_1", 2))
        self.assertTrue(buffered_data.HasValue("data"))

    def test_SetValue(self):
        buffered_data = BufferedDict(3)
        buffered_data.SetValue("data/sub_1", 1, 1)
        with self.assertRaises(RuntimeError):
            _ = buffered_data.GetValue("data/sub_1/int")

        with self.assertRaises(RuntimeError):
            _ = buffered_data.SetValue("data", 2)

    def test_GetMap(self):
        buffered_data = BufferedDict(3)
        buffered_data.SetValue("data/test/int", 1)
        buffered_data.SetValue("data/test/float", 3.0)
        buffered_data.SetValue("data/test/float", 4.0, 1)
        buffered_data.SetValue("data/test/str", "test", 1)
        buffered_data.SetValue("data/str", "old", 2)
        buffered_data.SetValue("data/int", 10)

        result = buffered_data.GetMap()
        for k, v in result.items():
            self.assertEqual(buffered_data[k], v)
        self.assertEqual(
            result,
            {
                "data": buffered_data.GetValue("data"),
                "data/int": 10,
                "data/test": buffered_data.GetValue("data/test"),
                "data/test/int": 1,
                "data/test/float": 3.0
            })

        result = buffered_data.GetMap(1)
        for k, v in result.items():
            self.assertEqual(buffered_data[k, 1], v)
        self.assertEqual(
            result,
            {
                "data": buffered_data.GetValue("data"),
                "data/test": buffered_data.GetValue("data/test"),
                "data/test/str": "test",
                "data/test/float": 4.0
            })

        result = buffered_data.GetMap(2)
        for k, v in result.items():
            self.assertEqual(buffered_data[k, 2], v)
        self.assertEqual(
            result,
            {
                "data": buffered_data.GetValue("data"),
                "data/str": "old",
                "data/test": buffered_data.GetValue("data/test")
            })

    def test_Buffers(self):
        buffered_data = BufferedDict(3)
        buffered_data.SetValue("data/sub_1/sub_sub", BufferedDict(4))
        buffered_data.SetValue("data/sub_2", BufferedDict(5))

        self.assertEqual(buffered_data.GetBufferSize(), 3)
        self.assertEqual(buffered_data.GetValue("data").GetBufferSize(), 3)
        self.assertEqual(buffered_data.GetValue("data/sub_1").GetBufferSize(), 3)
        self.assertEqual(buffered_data.GetValue("data/sub_1/sub_sub").GetBufferSize(), 4)
        self.assertEqual(buffered_data.GetValue("data/sub_2").GetBufferSize(), 5)

    def test_GetParent(self):
        buffered_data = BufferedDict(3)
        buffered_data.SetValue("data/sub_1/sub_sub", BufferedDict(4))
        buffered_data.SetValue("data/sub_2", BufferedDict(5))

        self.assertEqual(buffered_data.GetValue("data/sub_1/sub_sub").GetParent(), buffered_data.GetValue("data/sub_1"))
        self.assertEqual(buffered_data.GetValue("data/sub_1").GetParent(), buffered_data.GetValue("data"))
        self.assertEqual(buffered_data.GetValue("data/sub_2").GetParent(), buffered_data.GetValue("data"))

    def test_GetRoot(self):
        buffered_data = BufferedDict(3)
        buffered_data.SetValue("data/sub_1/sub_sub", BufferedDict(4))
        buffered_data.SetValue("data/sub_2", BufferedDict(5))

        self.assertEqual(buffered_data.GetValue("data/sub_1/sub_sub").GetRoot(), buffered_data)
        self.assertEqual(buffered_data.GetValue("data/sub_1").GetRoot(), buffered_data)
        self.assertEqual(buffered_data.GetValue("data/sub_2").GetRoot(), buffered_data)
        self.assertEqual(buffered_data.GetValue("data").GetRoot(), buffered_data)
        self.assertEqual(buffered_data.GetRoot(), buffered_data)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()