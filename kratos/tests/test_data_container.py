import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestDataContainer(KratosUnittest.TestCase):
    """Tests for the pybind11 bindings of Kratos::DataContainer."""

    def test_add_and_has(self):
        container = KM.DataContainer()

        self.assertFalse(container.Has(KM.PRESSURE))
        self.assertFalse(container.Has(KM.PRESSURE, KM.DoubleDataValuePolicy()))

        container.Add(KM.PRESSURE, KM.DoubleDataValuePolicy())

        self.assertTrue(container.Has(KM.PRESSURE))
        self.assertTrue(container.Has(KM.PRESSURE, KM.DoubleDataValuePolicy()))
        self.assertTrue(container.Has(KM.PRESSURE, KM.DoubleDataValuePolicy(), KM.NonHistoricalDataPolicy()))
        self.assertFalse(container.Has(KM.PRESSURE, KM.DoubleDataValuePolicy(), KM.HistoricalDataPolicy(KM.StepCategory.TimeStep, 3)))
        self.assertFalse(container.Has(KM.TEMPERATURE))

    def test_get_data_span_numpy_view(self):
        container = KM.DataContainer()
        accessor = container.Add(KM.PRESSURE, KM.DoubleDataValuePolicy())

        span = container.GetDataSpan(accessor)
        self.assertEqual(span.shape, (256,))  # default chunk size
        self.assertEqual(str(span.dtype), "float64")
        for value in span:
            self.assertEqual(value, 0.0)  # zero-initialized

        # Writing through the NumPy array modifies the stored data (view, not copy)
        span[:] = 1.5
        span_refetched = container.GetDataSpan(accessor)
        for value in span_refetched:
            self.assertEqual(value, 1.5)

    def test_chunk_size(self):
        container = KM.DataContainer(32)
        accessor = container.Add(KM.PRESSURE, KM.DoubleDataValuePolicy())
        self.assertEqual(len(container.GetDataSpan(accessor)), 32)

    def test_custom_zero_value(self):
        container = KM.DataContainer(16)
        accessor = container.Add(KM.STEP, KM.IntegerDataValuePolicy(-1))
        span = container.GetDataSpan(accessor)
        for value in span:
            self.assertEqual(value, -1)
        self.assertEqual(KM.IntegerDataValuePolicy(-1).Zero(), -1)

    def test_array3_span_shape(self):
        container = KM.DataContainer(8)
        # array_1d's default constructor does not zero-initialize its storage (unlike a
        # plain Python/NumPy sequence), so the policy zero must be passed explicitly here.
        accessor = container.Add(KM.VELOCITY, KM.Array1DDataValuePolicy3([0.0, 0.0, 0.0]))

        span = container.GetDataSpan(accessor)
        self.assertEqual(span.shape, (8, 3))
        self.assertTrue((span == 0.0).all())

        # Per-component write through the view
        span[:, 1] = 2.5
        span_refetched = container.GetDataSpan(accessor)
        self.assertTrue((span_refetched[:, 0] == 0.0).all())
        self.assertTrue((span_refetched[:, 1] == 2.5).all())
        self.assertTrue((span_refetched[:, 2] == 0.0).all())

    def test_string_span(self):
        container = KM.DataContainer(4)
        accessor = container.Add(KM.IDENTIFIER, KM.StringDataValuePolicy())

        span = container.GetDataSpan(accessor)
        self.assertEqual(len(span), 4)
        self.assertEqual(span[0], "")

        span[2] = "hello"
        self.assertEqual(span[2], "hello")
        self.assertEqual(span[-2], "hello")  # negative indices supported
        self.assertEqual(container.GetDataSpan(accessor).ToList(), ["", "", "hello", ""])

        with self.assertRaises(IndexError):
            span[4]

    def test_historical_steps(self):
        container = KM.DataContainer(8)
        history = KM.HistoricalDataPolicy(KM.StepCategory.TimeStep, 3)
        self.assertEqual(history.GetTotalNumberOfSteps(), 3)

        accessor = container.Add(KM.PRESSURE, KM.DoubleDataValuePolicy(), history)

        current = container.GetDataSpan(accessor)
        current[:] = 1.0

        previous_accessor = accessor.GetStepAccessor(KM.StepCategory.TimeStep, 1)
        previous = container.GetDataSpan(previous_accessor)
        for value in previous:
            self.assertEqual(value, 0.0)  # previous step untouched
        previous[:] = 2.0

        container.CloneStepData(KM.StepCategory.TimeStep)

        # After cloning, the new current step is a copy of the pre-clone current one
        self.assertTrue((container.GetDataSpan(accessor) == 1.0).all())
        self.assertTrue((container.GetDataSpan(previous_accessor) == 1.0).all())
        two_ago = container.GetDataSpan(accessor.GetStepAccessor(KM.StepCategory.TimeStep, 2))
        self.assertTrue((two_ago == 2.0).all())

    def test_accessor_equality(self):
        container = KM.DataContainer()
        accessor = container.Add(KM.PRESSURE, KM.DoubleDataValuePolicy())
        again = container.Add(KM.PRESSURE, KM.DoubleDataValuePolicy())

        self.assertTrue(accessor == again)
        self.assertEqual(accessor.GetIndex(), 0)
        self.assertFalse(accessor == accessor.GetStepAccessor(KM.StepCategory.TimeStep, 1))
        self.assertEqual(accessor.GetStepAccessor(KM.StepCategory.TimeStep, 1).GetStepBeforeCurrent(), 1)

    def test_sparse_workflow(self):
        # Mirror of the C++ sparse scenarios: 125 entities, activate the even ones (63),
        # then add the multiples of 5 (12 new entries -> 75)
        container = KM.DataContainer(125)
        index_accessor = container.Add(KM.STEP, KM.IntegerDataValuePolicy())
        value_accessor = container.Add(KM.PRESSURE, KM.DoubleSparseDataValuePolicy(index_accessor))

        self.assertTrue(KM.DoubleSparseDataValuePolicy(index_accessor).IsSparse())
        self.assertFalse(KM.DoubleDataValuePolicy().IsSparse())

        # Sparse chunks start empty
        self.assertEqual(len(container.GetDataSpan(value_accessor)), 0)

        index_span = container.GetDataSpan(index_accessor)
        counter = 0
        for i in range(len(index_span)):
            if i % 2 == 0:
                index_span[i] = counter
                counter += 1
            else:
                index_span[i] = -1

        container.UpdateSparseStorage(index_accessor)
        value_span = container.GetDataSpan(value_accessor)
        self.assertEqual(len(value_span), 125 // 2 + 1)

        value_span[:] = 1.0

        container.AddToSparseStorage(index_accessor, list(range(0, 125, 5)))

        value_span = container.GetDataSpan(value_accessor)  # re-fetch: storage was reallocated
        index_span = container.GetDataSpan(index_accessor)
        self.assertEqual(len(value_span), 125 // 2 + 13)

        for i in range(len(index_span)):
            if i % 2 == 0:
                self.assertEqual(value_span[index_span[i]], 1.0)  # original entries preserved
            elif i % 5 == 0:
                self.assertEqual(value_span[index_span[i]], 0.0)  # added entries zero-initialized
            else:
                self.assertEqual(index_span[i], -1)  # still inactive

    def test_resize(self):
        container = KM.DataContainer(4)
        accessor = container.Add(KM.PRESSURE, KM.DoubleDataValuePolicy())

        span = container.GetDataSpan(accessor)
        span[:] = [1.0, 2.0, 3.0, 4.0]

        container.Resize(6)

        span = container.GetDataSpan(accessor)  # re-fetch: the old view is invalidated
        self.assertEqual(len(span), 6)
        self.assertEqual(list(span[:4]), [1.0, 2.0, 3.0, 4.0])  # values preserved
        self.assertEqual(list(span[4:]), [0.0, 0.0])            # tail zero-initialized

        container.Resize(2)
        span = container.GetDataSpan(accessor)
        self.assertEqual(list(span), [1.0, 2.0])                # shrink keeps leading values

    def test_errors(self):
        container = KM.DataContainer()

        # Incompatible value policy: no matching Add overload exists for this combination
        with self.assertRaises(TypeError):
            container.Add(KM.STEP, KM.DoubleDataValuePolicy())

        # Missing data
        with self.assertRaises(RuntimeError):
            container.GetAccessor(KM.TEMPERATURE, KM.DoubleDataValuePolicy())

        # Existing variable with a different policy configuration
        container.Add(KM.PRESSURE, KM.DoubleDataValuePolicy())
        with self.assertRaises(RuntimeError):
            container.Add(KM.PRESSURE, KM.DoubleDataValuePolicy(5.0))


if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()
