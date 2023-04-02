
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestOptimizationInfo(kratos_unittest.TestCase):
    def test_SetGetValue(self):
        optimization_info = KratosOA.OptimizationInfo(1)

        optimization_info.SetValue("bool", True)
        optimization_info.SetValue("int", 1)
        optimization_info.SetValue("double", 2.0)
        optimization_info.SetValue("string", "test")

        sub_item = KratosOA.OptimizationInfo(3)
        sub_item.SetValue("sub_int", 3)
        sub_item.SetValue("sub_double", 4.0)
        optimization_info.SetValue("sub_item", sub_item)

        optimization_info.SetValue("sub_item2/sub_sub_item/int", 5)
        optimization_info.SetValue("sub_item2/double", 6.0)

        with self.assertRaises(RuntimeError):
            optimization_info.SetValue("sub_item2/double/double", 7.0)

        self.assertEqual(optimization_info.GetBool("bool"), True)

        with self.assertRaises(RuntimeError):
            optimization_info.GetString("bool")

        self.assertEqual(optimization_info.GetInt("int"), 1)
        self.assertEqual(optimization_info.GetDouble("double"), 2.0)
        self.assertEqual(optimization_info.GetString("string"), "test")
        self.assertEqual(optimization_info.GetInt("sub_item/sub_int"), 3)
        self.assertEqual(optimization_info.GetDouble("sub_item2/double"), 6)
        self.assertEqual(optimization_info.GetInt("sub_item2/sub_sub_item/int"), 5)

    def test_SetBufferSize(self):
        optimization_info = KratosOA.OptimizationInfo(1)
        optimization_info.SetValue("int", 1)
        optimization_info.SetValue("sub_item/int", 2)

        sub_item = optimization_info.GetSubItem("sub_item")
        self.assertEqual(sub_item.GetBufferSize(), 1)

        optimization_info.SetBufferSize(2)
        self.assertEqual(optimization_info.GetBufferSize(), 2)
        self.assertEqual(sub_item.GetBufferSize(), 1)

        optimization_info.SetBufferSize(2, True)
        self.assertEqual(optimization_info.GetBufferSize(), 2)
        self.assertEqual(sub_item.GetBufferSize(), 2)

    def test_IsValue(self):
        optimization_info = KratosOA.OptimizationInfo(1)

        optimization_info.SetValue("bool", True)
        optimization_info.SetValue("int", 1)
        optimization_info.SetValue("double", 2.0)
        optimization_info.SetValue("string", "test")

        self.assertTrue(optimization_info.IsBool("bool"))
        self.assertTrue(optimization_info.IsInt("int"))
        self.assertTrue(optimization_info.IsDouble("double"))
        self.assertTrue(optimization_info.IsString("string"))

    def test_Buffer(self):
        optimization_info = KratosOA.OptimizationInfo(2)

        optimization_info.SetValue("bool", True)
        optimization_info.SetValue("string", "test")
        optimization_info.SetValue("string", "old", 1)

        sub_item = KratosOA.OptimizationInfo(3)
        optimization_info.SetValue("sub_item", sub_item)
        optimization_info.SetValue("sub_item/int", 1)
        optimization_info.SetValue("sub_item/int", 2, 1)
        optimization_info.SetValue("sub_item/int", 3, 2)

        self.assertEqual(optimization_info.GetBool("bool", 0), True)
        self.assertEqual(optimization_info.GetString("string", 0), "test")
        self.assertEqual(optimization_info.GetInt("sub_item/int", 0), 1)
        self.assertEqual(optimization_info.GetInt("sub_item/int", 1), 2)
        self.assertEqual(optimization_info.GetInt("sub_item/int", 2), 3)

        optimization_info.SetValue("bool", "old", 1)

        self.assertTrue(optimization_info.IsBool("bool", 0))
        self.assertTrue(optimization_info.IsString("bool", 1))

        sub_item.SetValue("hello/hello/hello", "hello", 1)
        sub_item.SetValue("hello/hello/step", "step", 2)
        self.assertTrue(sub_item.HasValue("hello"))
        self.assertTrue(sub_item.HasValue("hello/hello"))
        self.assertFalse(sub_item.HasValue("hello/hello/hello"))
        self.assertTrue(sub_item.HasValue("hello/hello/hello", 1))
        self.assertFalse(sub_item.HasValue("hello/hello/step", 1))
        self.assertTrue(sub_item.HasValue("hello/hello/step", 2))

    def test_AdvanceStep(self):
        optimization_info = KratosOA.OptimizationInfo(2)

        optimization_info.SetValue("bool", True)
        optimization_info.SetValue("bool", "old", 1)
        optimization_info.SetValue("string", "test")
        optimization_info.SetValue("string", "old", 1)

        sub_item = KratosOA.OptimizationInfo(3)
        optimization_info.SetValue("sub_item", sub_item)
        optimization_info.SetValue("sub_item/int", 1)
        optimization_info.SetValue("sub_item/int", 2, 1)
        optimization_info.SetValue("sub_item/int", 3, 2)

        optimization_info.AdvanceStep()
        self.assertEqual(optimization_info.GetString("bool", 0), "old")
        self.assertEqual(optimization_info.GetBool("bool", 1), True)
        self.assertEqual(optimization_info.GetString("string", 0), "old")
        self.assertEqual(optimization_info.GetString("string", 1), "test")
        self.assertEqual(optimization_info.GetInt("sub_item/int", 0), 2)
        self.assertEqual(optimization_info.GetInt("sub_item/int", 1), 3)
        self.assertEqual(optimization_info.GetInt("sub_item/int", 2), 1)

        optimization_info.AdvanceStep()
        self.assertEqual(optimization_info.GetString("bool", 1), "old")
        self.assertEqual(optimization_info.GetBool("bool", 0), True)
        self.assertEqual(optimization_info.GetString("string", 1), "old")
        self.assertEqual(optimization_info.GetString("string", 0), "test")
        self.assertEqual(optimization_info.GetInt("sub_item/int", 0), 3)
        self.assertEqual(optimization_info.GetInt("sub_item/int", 1), 1)
        self.assertEqual(optimization_info.GetInt("sub_item/int", 2), 2)

    def test_Misc(self):
        optimization_info = KratosOA.OptimizationInfo(2)
        optimization_info.SetValue("test/test1/int", 1)

        with self.assertRaises(RuntimeError):
            optimization_info.SetValue("test/test1", 10)

        with self.assertRaises(RuntimeError):
            optimization_info.SetValue("test/test1/int/hello", "hello")

        self.assertFalse(optimization_info.IsInt("test/test1"))
        optimization_info.SetValue("test/test1", 10, overwrite=True)
        self.assertTrue(optimization_info.IsInt("test/test1"))

        self.assertFalse(optimization_info.IsSubItem("test/test1"))
        optimization_info.SetValue("test/test1", KratosOA.OptimizationInfo(3), overwrite=True)
        self.assertTrue(optimization_info.IsSubItem("test/test1"))

    def test_SetContainers(self):
        optimization_info = KratosOA.OptimizationInfo(2)
        optimization_info.SetValue("test/test1/int", 1)

        model = Kratos.Model()
        model_part = model.CreateModelPart("test")

        a = Kratos.ContainerExpression.HistoricalExpression(model_part)
        optimization_info.SetValue("test/container", a)
        self.assertTrue(optimization_info.IsHistoricalExpression("test/container"))
        self.assertEqual(optimization_info.GetHistoricalExpression("test/container"), a)

        optimization_info.SetValue("test/model_part", model_part)
        self.assertTrue(optimization_info.IsModelPart("test/model_part"))

        c = optimization_info.GetModelPart("test/model_part")
        self.assertEqual(c, model_part)

    def test_GetKeys(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        a = Kratos.ContainerExpression.HistoricalExpression(model_part)

        optimization_info = KratosOA.OptimizationInfo(2)

        optimization_info.SetValue("bool", True)
        optimization_info.SetValue("int", 1)
        optimization_info.SetValue("double", 2.0)
        optimization_info.SetValue("string", "test")

        sub_item = KratosOA.OptimizationInfo(3)
        sub_item.SetValue("sub_int", 3)
        sub_item.SetValue("sub_double", 4.0)
        optimization_info.SetValue("sub_item", sub_item)

        optimization_info.SetValue("sub_item2/sub_sub_item/int", 5)
        optimization_info.SetValue("sub_item2/double", 6.0)

        sub_item.SetValue("sub_int", 3, 1)
        sub_item.SetValue("sub_double", 4.0, 1)

        optimization_info.SetValue("sub_item/sub_model_part/model_part", model_part, 2)
        optimization_info.SetValue("sub_item/container/contaner", a, 2)

        for k in optimization_info.GetKeys():
            self.assertTrue(optimization_info.HasValue(k))

        for k in optimization_info.GetKeys(1):
            self.assertTrue(optimization_info.HasValue(k, 1))

        for k in optimization_info.GetKeys(2):
            self.assertTrue(optimization_info.HasValue(k, 2))

        self.assertEqual(sorted(optimization_info.GetKeys()), sorted(["bool", "int", "double", "string", "sub_item/sub_int", "sub_item/sub_double", "sub_item2/sub_sub_item/int", "sub_item2/double"]))
        self.assertEqual(sorted(optimization_info.GetKeys(1)), sorted(["sub_item/sub_int", "sub_item/sub_double"]))
        self.assertEqual(sorted(optimization_info.GetKeys(2)), sorted(["sub_item/sub_model_part/model_part", "sub_item/container/contaner"]))

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.TESTS_OUTPUTS)  # TESTS_OUTPUTS
    kratos_unittest.main()