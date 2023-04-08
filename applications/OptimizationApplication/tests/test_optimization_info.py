
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo

class TestOptimizationInfo(kratos_unittest.TestCase):
    def test_AdvanceStep(self):
        optimization_info = OptimizationInfo(3)

        optimization_info["step"] = 1
        self.assertEqual(optimization_info["step"], 1)

        optimization_info.AdvanceStep()
        optimization_info.ClearStep()
        optimization_info["step"] = 2
        self.assertEqual(optimization_info["step"], 2)
        self.assertEqual(optimization_info["step", 1], 1)

        optimization_info.AdvanceStep()
        optimization_info.ClearStep()
        optimization_info["step"] = 3
        self.assertEqual(optimization_info["step"], 3)
        self.assertEqual(optimization_info["step", 1], 2)
        self.assertEqual(optimization_info["step", 2], 1)

        optimization_info.AdvanceStep()
        optimization_info.ClearStep()
        optimization_info["step"] = 4
        self.assertEqual(optimization_info["step"], 4)
        self.assertEqual(optimization_info["step", 1], 3)
        self.assertEqual(optimization_info["step", 2], 2)
        with self.assertRaises(RuntimeError):
            optimization_info["step", 3]

    def test_AdvanceStepSet(self):
        optimization_info = OptimizationInfo(3)
        optimization_info["test_1/test_sub_1"] = OptimizationInfo(2)

        self.assertEqual(optimization_info["test_1"].GetBufferSize(), 3)
        self.assertEqual(optimization_info["test_1/test_sub_1"].GetBufferSize(), 2)

        optimization_info["test_1/int"] = 1
        optimization_info["test_1/test_sub_1/int"] = 2

        optimization_info.AdvanceStep()
        optimization_info.ClearStep()
        optimization_info["test_1/int"] = 3
        optimization_info["test_1/test_sub_1/int"] = 4

        optimization_info.AdvanceStep()
        optimization_info.ClearStep()
        optimization_info["test_1/int"] = 5
        optimization_info["test_1/test_sub_1/int"] = 6

        optimization_info.AdvanceStep()
        optimization_info.ClearStep()
        optimization_info["test_1/int"] = 7
        optimization_info["test_1/test_sub_1/int"] = 8

        optimization_info.AdvanceStep()
        optimization_info.ClearStep()
        optimization_info["test_1/int"] = 9
        optimization_info["test_1/test_sub_1/int"] = 10

        optimization_info.AdvanceStep()
        optimization_info.ClearStep()
        optimization_info["test_1/int"] = 11
        optimization_info["test_1/test_sub_1/int"] = 12

        self.assertEqual(optimization_info["test_1/int", 0], 11)
        self.assertEqual(optimization_info["test_1/int", 1], 9)
        self.assertEqual(optimization_info["test_1/int", 2], 7)
        with self.assertRaises(RuntimeError):
            optimization_info["test_1/int", 3]

        self.assertEqual(optimization_info["test_1/test_sub_1/int", 0], 12)
        self.assertEqual(optimization_info["test_1/test_sub_1/int", 1], 10)
        with self.assertRaises(RuntimeError):
            optimization_info["test_1/test_sub_1/int", 3]

    def test_HasValue(self):
        optimization_info = OptimizationInfo(2)
        optimization_info["data/sub_1"] = 1
        optimization_info["data/sub_2", 1] = 2
        optimization_info.SetValue("data/sub_3/sub_sub1", 3, 1)

        self.assertTrue(optimization_info.HasValue("data/sub_2", 1))
        self.assertTrue(optimization_info["data"].HasValue("sub_3/sub_sub1", 1))
        self.assertTrue(optimization_info.HasValue("data/sub_2", 1))
        self.assertTrue(optimization_info.HasValue("data/sub_3", 0))
        self.assertTrue(optimization_info.HasValue("data/sub_3", 1))
        self.assertFalse(optimization_info.HasValue("data/sub_2"))
        self.assertFalse(optimization_info.HasValue("data_2/sub_2"))

    def test_RemoveValue(self):
        optimization_info = OptimizationInfo(3)
        optimization_info["data/sub_1"] = 1
        optimization_info["data/sub_1", 1] = 2
        optimization_info["data/sub_1", 2] = 3

        self.assertTrue(optimization_info.HasValue("data/sub_1", 0))
        self.assertTrue(optimization_info.HasValue("data/sub_1", 1))
        self.assertTrue(optimization_info.HasValue("data/sub_1", 2))

        del optimization_info["data/sub_1", 0]
        self.assertFalse(optimization_info.HasValue("data/sub_1", 0))
        self.assertTrue(optimization_info.HasValue("data/sub_1", 1))
        self.assertTrue(optimization_info.HasValue("data/sub_1", 2))

        del optimization_info["data/sub_1", 1]
        self.assertFalse(optimization_info.HasValue("data/sub_1", 0))
        self.assertFalse(optimization_info.HasValue("data/sub_1", 1))
        self.assertTrue(optimization_info.HasValue("data/sub_1", 2))

        del optimization_info["data/sub_1", 2]
        self.assertFalse(optimization_info.HasValue("data/sub_1", 0))
        self.assertFalse(optimization_info.HasValue("data/sub_1", 1))
        self.assertFalse(optimization_info.HasValue("data/sub_1", 2))
        self.assertTrue(optimization_info.HasValue("data"))

    def test_SetValue(self):
        optimization_info = OptimizationInfo(3)
        optimization_info["data/sub_1", 1] = 1
        with self.assertRaises(RuntimeError):
            optimization_info["data/sub_1/int"]

        with self.assertRaises(RuntimeError):
            optimization_info["data"] = 2

    def test_GetMap(self):
        optimization_info = OptimizationInfo(3)
        optimization_info["data/test/int"] = 1
        optimization_info["data/test/float"] = 3.0
        optimization_info["data/test/float", 1] = 4.0
        optimization_info["data/test/str", 1] = "test"
        optimization_info["data/str", 2] = "old"
        optimization_info["data/int"] = 10

        result = optimization_info.GetMap()
        for k, v in result.items():
            self.assertEqual(optimization_info[k], v)
        self.assertEqual(
            result,
            {
                "data": optimization_info["data"],
                "data/int": 10,
                "data/test": optimization_info["data/test"],
                "data/test/int": 1,
                "data/test/float": 3.0
            })

        result = optimization_info.GetMap(1)
        for k, v in result.items():
            self.assertEqual(optimization_info[k, 1], v)
        self.assertEqual(
            result,
            {
                "data": optimization_info["data"],
                "data/test": optimization_info["data/test"],
                "data/test/str": "test",
                "data/test/float": 4.0
            })

        result = optimization_info.GetMap(2)
        for k, v in result.items():
            self.assertEqual(optimization_info[k, 2], v)
        self.assertEqual(
            result,
            {
                "data": optimization_info["data"],
                "data/str": "old",
                "data/test": optimization_info["data/test"]
            })

    def test_Buffers(self):
        optimization_info = OptimizationInfo(3)
        optimization_info["data/sub_1/sub_sub"] = OptimizationInfo(4)
        optimization_info["data/sub_2"] = OptimizationInfo(5)

        self.assertEqual(optimization_info.GetBufferSize(), 3)
        self.assertEqual(optimization_info["data"].GetBufferSize(), 3)
        self.assertEqual(optimization_info["data/sub_1"].GetBufferSize(), 3)
        self.assertEqual(optimization_info["data/sub_1/sub_sub"].GetBufferSize(), 4)
        self.assertEqual(optimization_info["data/sub_2"].GetBufferSize(), 5)

    def test_GetParent(self):
        optimization_info = OptimizationInfo(3)
        optimization_info["data/sub_1/sub_sub"] = OptimizationInfo(4)
        optimization_info["data/sub_2"] = OptimizationInfo(5)

        self.assertEqual(optimization_info["data/sub_1/sub_sub"].GetParent(), optimization_info["data/sub_1"])
        self.assertEqual(optimization_info["data/sub_1"].GetParent(), optimization_info["data"])
        self.assertEqual(optimization_info["data/sub_2"].GetParent(), optimization_info["data"])

    def test_GetRoot(self):
        optimization_info = OptimizationInfo(3)
        optimization_info["data/sub_1/sub_sub"] = OptimizationInfo(4)
        optimization_info["data/sub_2"] = OptimizationInfo(5)

        self.assertEqual(optimization_info["data/sub_1/sub_sub"].GetRoot(), optimization_info)
        self.assertEqual(optimization_info["data/sub_1"].GetRoot(), optimization_info)
        self.assertEqual(optimization_info["data/sub_2"].GetRoot(), optimization_info)
        self.assertEqual(optimization_info["data"].GetRoot(), optimization_info)
        self.assertEqual(optimization_info.GetRoot(), optimization_info)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()