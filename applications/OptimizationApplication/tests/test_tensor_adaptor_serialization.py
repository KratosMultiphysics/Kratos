import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.OptimizationApplication as KratosOA

class TestPropertiesVariableTensorAdaptorSerialization(kratos_unittest.TestCase):
    def test_CombinedTensorAdaptorWithPropertiesVariableChildUniqueProperties(self):
        model = Kratos.Model()
        model_part = self._CreateTriangleElements(model, 2, shared_properties=False)

        original = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(model_part.Elements, Kratos.DENSITY)
        original.CollectData()
        original.ViewData()[:] = [11.0, 22.0]
        original.StoreData()
        original.Check()

        loaded_model_part = self._CheckCombinedPropertiesTensorAdaptorRoundTrip(model, original)
        self.assertEqual(loaded_model_part.NumberOfProperties(), 2)
        self.assertEqual([element.Properties.GetValue(Kratos.DENSITY) for element in loaded_model_part.Elements], [11.0, 22.0])

    def test_CombinedTensorAdaptorWithPropertiesVariableChildSharedProperties(self):
        # many elements on one Properties without DENSITY: StoreData() must create it once, not race on it
        model = Kratos.Model()
        model_part = self._CreateTriangleElements(model, 100, shared_properties=True)
        properties = model_part.GetProperties(1)
        self.assertFalse(properties.Has(Kratos.DENSITY))

        original = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(model_part.Elements, Kratos.DENSITY)
        original.ViewData()[:] = 11.0
        original.StoreData()
        self.assertEqual(properties.GetValue(Kratos.DENSITY), 11.0)

        # entities sharing a Properties cannot hold different values
        original.ViewData()[57] = 12.0
        with self.assertRaisesRegex(RuntimeError, "DENSITY"):
            original.StoreData()
        self.assertEqual(properties.GetValue(Kratos.DENSITY), 11.0)
        original.ViewData()[57] = 11.0

        loaded_model_part = self._CheckCombinedPropertiesTensorAdaptorRoundTrip(model, original)
        self.assertEqual(loaded_model_part.NumberOfProperties(), 1)
        self.assertTrue(all(element.Properties.Id == 1 for element in loaded_model_part.Elements))
        self.assertEqual(loaded_model_part.GetProperties(1).GetValue(Kratos.DENSITY), 11.0)

    @staticmethod
    def _CreateTriangleElements(model: Kratos.Model, number_of_elements: int, shared_properties: bool) -> Kratos.ModelPart:
        model_part = model.CreateModelPart("Test")
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        properties = model_part.CreateNewProperties(1)
        for i in range(number_of_elements):
            if not shared_properties and i > 0:
                properties = model_part.CreateNewProperties(i + 1)
            model_part.CreateNewElement("Element2D3N", i + 1, [1, 2, 3], properties)
        return model_part

    def _CheckCombinedPropertiesTensorAdaptorRoundTrip(self, model: Kratos.Model, original: KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor) -> Kratos.ModelPart:
        # children are saved through TensorAdaptor<double> pointers, so PropertiesVariableTensorAdaptor
        # must be registered with the Serializer (optimization_application.cpp)
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

        self._CheckValues(loaded.data, combined.data)
        loaded_children = loaded.GetTensorAdaptors()
        self.assertEqual(len(loaded_children), 1)
        self.assertIsInstance(loaded_children[0], KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor)

        # the loaded adaptor writes back into the loaded model part
        loaded_children[0].StoreData()
        return load_model.GetModelPart("Test")

    def test_CombinedTensorAdaptorWithVelocityChildren(self):
        for copy in [True, False]:
            with self.subTest(copy=copy):
                self._CheckCombinedVelocityTensorAdaptorRoundTrip(copy)

    def _CheckCombinedVelocityTensorAdaptorRoundTrip(self, copy: bool) -> None:
        model = Kratos.Model()
        model_part = model.CreateModelPart("Test")
        for i in range(3):
            model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0).SetValue(Kratos.VELOCITY, Kratos.Array3([i + 1.0, i + 2.0, i + 3.0]))
        for i in range(2):
            # one Properties per element, so StoreData() writes independent values
            properties = model_part.CreateNewProperties(i + 1)
            properties.SetValue(Kratos.VELOCITY, Kratos.Array3([10.0 * (i + 1), 20.0 * (i + 1), 30.0 * (i + 1)]))
            model_part.CreateNewElement("Element2D3N", i + 1, [1, 2, 3], properties)

        properties_ta = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(model_part.Elements, Kratos.VELOCITY)
        nodal_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(model_part.Nodes, Kratos.VELOCITY)
        combined = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([properties_ta, nodal_ta], 0, perform_collect_data_recursively=True, perform_store_data_recursively=True, copy=copy)
        combined.CollectData()
        original_data = numpy.array(combined.data)
        self.assertEqual(list(combined.Shape()), [5, 3])
        self._CheckValues(original_data, [[10, 20, 30], [20, 40, 60], [1, 2, 3], [2, 3, 4], [3, 4, 5]])

        # the copy flag itself is not saved, only its effect: copy=False children share the storage of
        # the given adaptors. Saving those too lets the Serializer's pointer tracking restore the sharing.
        serializer = Kratos.StreamSerializer()
        serializer.Set(Kratos.Serializer.SHALLOW_GLOBAL_POINTERS_SERIALIZATION)
        serializer.Save("Model", model)
        serializer.Save("Combined", combined)
        serializer.Save("Properties", properties_ta)
        serializer.Save("Nodal", nodal_ta)

        load_model = Kratos.Model()
        loaded = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor()
        loaded_properties_ta = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor()
        loaded_nodal_ta = Kratos.TensorAdaptors.VariableTensorAdaptor()
        serializer.Load("Model", load_model)
        serializer.Load("Combined", loaded)
        serializer.Load("Properties", loaded_properties_ta)
        serializer.Load("Nodal", loaded_nodal_ta)

        self.assertEqual(list(loaded.Shape()), [5, 3])
        self._CheckValues(loaded.data, original_data)
        loaded_children = loaded.GetTensorAdaptors()
        self.assertEqual(len(loaded_children), 2)
        self.assertIsInstance(loaded_children[0], KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor)
        self.assertIsInstance(loaded_children[1], Kratos.TensorAdaptors.VariableTensorAdaptor)
        self.assertEqual(list(loaded_children[0].Shape()), [2, 3])
        self.assertEqual(list(loaded_children[1].Shape()), [3, 3])

        original_children = combined.GetTensorAdaptors()
        for original_child, loaded_child, given, loaded_given in zip(original_children, loaded_children, [properties_ta, nodal_ta], [loaded_properties_ta, loaded_nodal_ta]):
            self.assertEqual(numpy.shares_memory(original_child.data, given.data), not copy)
            self.assertEqual(numpy.shares_memory(loaded_child.data, loaded_given.data), not copy)

        # recursive StoreData() reaches the loaded model's entities, not the original model's
        loaded.data[:] = -original_data
        loaded.StoreData()
        load_model_part = load_model.GetModelPart("Test")
        self._CheckValues([element.Properties[Kratos.VELOCITY] for element in load_model_part.Elements], -original_data[:2])
        self._CheckValues([node.GetValue(Kratos.VELOCITY) for node in load_model_part.Nodes], -original_data[2:])
        self._CheckValues([element.Properties[Kratos.VELOCITY] for element in model_part.Elements], original_data[:2])
        self._CheckValues([node.GetValue(Kratos.VELOCITY) for node in model_part.Nodes], original_data[2:])

        # recursive CollectData() reads the loaded model's entities into the children and the combined data
        for node in load_model_part.Nodes:
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([0.5 * node.Id, 0.0, 0.0]))
        loaded.CollectData()
        self._CheckValues(loaded_children[1].data, [[0.5, 0, 0], [1.0, 0, 0], [1.5, 0, 0]])
        self._CheckValues(loaded.data[2:], [[0.5, 0, 0], [1.0, 0, 0], [1.5, 0, 0]])
        self._CheckValues(loaded.data[:2], -original_data[:2])
        self._CheckValues(combined.data, original_data)

    def _CheckValues(self, values, expected) -> None:
        self.assertVectorAlmostEqual(numpy.ravel(numpy.array(values)), numpy.ravel(numpy.array(expected, dtype=float)))

if __name__ == "__main__":
    kratos_unittest.main()
