import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.OptimizationApplication as KratosOA

class TestPropertiesVariableTensorAdaptorSerialization(kratos_unittest.TestCase):
    def test_CombinedTensorAdaptorWithPropertiesVariableChild(self):
        # children are saved through TensorAdaptor<double> pointers, so PropertiesVariableTensorAdaptor
        # must be registered with the Serializer (optimization_application.cpp), otherwise saving
        # throws "no object registered with type id ...".
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
