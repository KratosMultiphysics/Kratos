import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.OptimizationApplication as KratosOA

class TestPropertiesVariableTensorAdaptorSerialization(kratos_unittest.TestCase):
    def test_CombinedTensorAdaptorWithPropertiesVariableChild(self):
        # Regression: PropertiesVariableTensorAdaptor had no default constructor, serialization
        # hooks, or Serializer::Register entry, so saving a combined field whose child survived
        # Clone() as this (real, OptimizationApplication) type threw "no object registered with
        # type id ...", the same class of bug this PR's core subtype registrations fixed.
        model = Kratos.Model()
        model_part = model.CreateModelPart("Test")
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        properties = model_part.CreateNewProperties(1)
        model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)
        model_part.CreateNewElement("Element2D3N", 2, [1, 2, 3], properties)

        original = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(model_part.Elements, Kratos.DENSITY)
        original.CollectData()
        original.ViewData()[:] = [11.0, 22.0]
        original.StoreData()

        combined = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([original], False, False, False)
        combined.CollectData()

        serializer = Kratos.StreamSerializer()
        serializer.Set(Kratos.Serializer.SHALLOW_GLOBAL_POINTERS_SERIALIZATION)
        serializer.Save("Model", model)
        serializer.Save("TA", combined)

        load_model = Kratos.Model()
        loaded = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor()
        serializer.Load("Model", load_model)
        serializer.Load("TA", loaded)

        self.assertEqual(loaded.Size(), combined.Size())
        for i in range(combined.Size()):
            self.assertAlmostEqual(loaded.ViewData()[i], combined.ViewData()[i])

        loaded_children = loaded.GetTensorAdaptors()
        self.assertEqual(len(loaded_children), 1)
        self.assertIsInstance(loaded_children[0], KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor)

if __name__ == "__main__":
    kratos_unittest.main()
