
import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as KratosUtilities
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestVariableTensorAdaptors(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.SetBufferSize(2)

        cls.input_list_of_variables = [
                Kratos.PRESSURE,                # double
                Kratos.VELOCITY,                # array3
                Kratos.EXTERNAL_FORCES_VECTOR,  # vector
                Kratos.NORMAL_SHAPE_DERIVATIVE  # matrix
            ]

        cls.output_list_of_variables = [
                Kratos.DENSITY,                 # double
                Kratos.ACCELERATION,            # array3
                Kratos.INTERNAL_FORCES_VECTOR,  # vector
                Kratos.PK2_STRESS_TENSOR        # matrix
            ]

        for var in cls.input_list_of_variables:
            cls.model_part.AddNodalSolutionStepVariable(var)

        for var in cls.output_list_of_variables:
            cls.model_part.AddNodalSolutionStepVariable(var)

        for i in range(10):
            cls.model_part.CreateNewNode(i + 1, i, i + 1, i + 2)

        for i in range(20):
            cls.model_part.CreateNewProperties(i + 1)

        for i in range(10):
            cls.model_part.CreateNewCondition("LineCondition2D2N", i + 1, [i + 1, ((i + 1) % 10 + 1)], cls.model_part.GetProperties(i + 1))

        for i in range(10):
            cls.model_part.CreateNewElement("Element3D2N", i + 1, [i + 1, ((i + 1) % 10 + 1)], cls.model_part.GetProperties(i + 1))

    def setUp(self):
        Kratos.VariableUtils().ClearNonHistoricalData(self.model_part.Nodes)
        Kratos.VariableUtils().ClearNonHistoricalData(self.model_part.Conditions)
        Kratos.VariableUtils().ClearNonHistoricalData(self.model_part.Elements)

        # # clear the values at the output vars
        # for var in self.output_list_of_variables:
        #     Kratos.VariableUtils().SetHistoricalVariableToZero(var, self.model_part.Nodes)

        # # clear the values at the input vars
        # for var in self.input_list_of_variables:
        #     Kratos.VariableUtils().SetHistoricalVariableToZero(var, self.model_part.Nodes)

        def setter(container, setter_method):
            for entity in container:
                setter_method(entity, Kratos.PRESSURE, entity.Id + 1)
                setter_method(entity, Kratos.VELOCITY, Kratos.Vector([entity.Id * 2 + 1, entity.Id * 3 + 1, entity.Id * 4 + 1]))
                setter_method(entity, Kratos.EXTERNAL_FORCES_VECTOR, Kratos.Vector([entity.Id * 2 + 1, entity.Id * 3 + 1, entity.Id * 4 + 1, entity.Id * 2 + 1, entity.Id * 3 + 1]))
                setter_method(entity, Kratos.NORMAL_SHAPE_DERIVATIVE, Kratos.Matrix([[entity.Id * 2 + 1, entity.Id * 3 + 1], [entity.Id * 4 + 1, entity.Id * 2 + 1], [entity.Id * 3 + 1, entity.Id * 4 + 1]]))

        setter(self.model_part.Nodes, lambda x, y, z: x.SetSolutionStepValue(y, 0, z * 1.5))
        setter(self.model_part.Nodes, lambda x, y, z: x.SetSolutionStepValue(y, 1, z * 1.5))
        setter(self.model_part.Nodes, lambda x, y, z: x.SetValue(y, z * 2.1))
        setter(self.model_part.Properties, lambda x, y, z: x.SetValue(y, z * 6.1))
        setter(self.model_part.Conditions, lambda x, y, z: x.SetValue(y, z * 3.1))
        setter(self.model_part.Elements, lambda x, y, z: x.SetValue(y, z * 4.1))

    def test_Shape1(self):
        t_adaptor_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.VELOCITY)
        self.assertVectorAlmostEqual(t_adaptor_1.Shape(), [10, 3])
        self.assertVectorAlmostEqual(t_adaptor_1.GetDataShape(), [3])

        t_adaptor_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.NORMAL_SHAPE_DERIVATIVE)
        self.assertVectorAlmostEqual(t_adaptor_1.Shape(), [10, 3, 2])
        self.assertVectorAlmostEqual(t_adaptor_1.GetDataShape(), [3, 2])

        t_adaptor_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.PK2_STRESS_TENSOR)
        self.assertVectorAlmostEqual(t_adaptor_1.Shape(), [10, 0, 0])
        self.assertVectorAlmostEqual(t_adaptor_1.GetDataShape(), [0, 0])

        t_adaptor_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.PK2_STRESS_TENSOR, data_shape=[5,6])
        self.assertVectorAlmostEqual(t_adaptor_1.Shape(), [10, 5, 6])
        self.assertVectorAlmostEqual(t_adaptor_1.GetDataShape(), [5, 6])

    def test_Shape2(self):
        t_adaptor_1 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)
        t_adaptor_1.CollectData()

        t_adaptor_2 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY, data_shape=[2])
        t_adaptor_2.CollectData()

        t_adaptor_3 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY, data_shape=[1])
        t_adaptor_3.CollectData()

        t_adaptor_4 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY_X)
        t_adaptor_4.CollectData()

        self.assertVectorAlmostEqual(t_adaptor_4.data, t_adaptor_3.data)
        self.assertVectorAlmostEqual(t_adaptor_3.data, t_adaptor_1.data[:, 0])
        for v1, v2 in zip(t_adaptor_2.data, t_adaptor_1.data[:, :2]):
            self.assertVectorAlmostEqual(v1, v2)

    def test_ScopedView(self):
        t_adaptor_1 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY_X)
        t_adaptor_1.CollectData()

        def temp(tensor_adaptor: Kratos.TensorAdaptors.DoubleTensorAdaptor):
            temp_numpy = tensor_adaptor.data
            new_numpy = temp_numpy * 2
            return new_numpy[:]

        test = temp(t_adaptor_1)
        for v1, v2 in zip(test, t_adaptor_1.data):
            self.assertEqual(v1, v2 * 2)

    def test_Assignment(self):
        t_adaptor_1 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY_X)
        t_adaptor_1.CollectData()

        t_adaptor_2 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.DENSITY)
        t_adaptor_2.CollectData()

        t_adaptor_2.data = t_adaptor_1.data
        t_adaptor_2.StoreData()
        numpy_orig = t_adaptor_2.data

        for node in self.model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(Kratos.VELOCITY_X), node.GetSolutionStepValue(Kratos.DENSITY))

        t_adaptor_3 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY_Y)
        t_adaptor_3.data = numpy_orig * 3
        t_adaptor_3.StoreData()

        for node in self.model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(Kratos.VELOCITY_Y), node.GetSolutionStepValue(Kratos.DENSITY) * 3)

    def test_NonInitializedVars(self):
        t_adaptor_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.INITIAL_STRAIN)
        numpy_data = t_adaptor_1.data

        t_adaptor_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.RECOVERED_STRESS)
        t_adaptor_2.data = t_adaptor_1.data

        t_adaptor_1.CollectData()
        t_adaptor_2.data = t_adaptor_1.data
        t_adaptor_2.StoreData()

    def test_SupportedNdArrays(self):
        t_adaptor = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.DENSITY)

        for i, numpy_dtype in enumerate([numpy.bool, numpy.uint8, numpy.uint16, numpy.uint32, numpy.uint64, numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64, numpy.float128]):
            numpy_array = numpy.array(numpy.random.random((t_adaptor.Shape())), dtype=numpy_dtype)
            t_adaptor.data = numpy_array
            t_adaptor.StoreData()

            # now check whether values are written properly
            for i, node in enumerate(self.model_part.Nodes):
                self.assertAlmostEqual(node.GetValue(Kratos.DENSITY), numpy_array[i])

        # non supported
        with self.assertRaises(RuntimeError):
            numpy_array = numpy.zeros(shape=(t_adaptor.Shape()), dtype=numpy.float16)
            t_adaptor.data = numpy_array

    def test_NodeVariableTensorAdaptor(self):
        self.__TestVariableTensorAdaptor(self.model_part.Nodes)

    def test_ConditionVariableTensorAdaptor(self):
        self.__TestVariableTensorAdaptor(self.model_part.Conditions)

    def test_ElementVariableTensorAdaptor(self):
        self.__TestVariableTensorAdaptor(self.model_part.Elements)

    def test_PropertyVariableTensorAdaptor(self):
        self.__TestVariableTensorAdaptor(self.model_part.Properties)

    def __TestVariableTensorAdaptor(self, container):
        for input_var, output_var in zip(self.input_list_of_variables, self.output_list_of_variables):
            read_tensor_adaptor = Kratos.TensorAdaptors.VariableTensorAdaptor(container, input_var)
            read_tensor_adaptor.CollectData()

            write_tensor_adaptor = Kratos.TensorAdaptors.VariableTensorAdaptor(read_tensor_adaptor.GetContainer(), output_var, data_shape=read_tensor_adaptor.GetDataShape())

            # modify the write tensor data
            write_tensor_adaptor.data = read_tensor_adaptor.data * 2

            write_tensor_adaptor.StoreData()
            for entity in container:
                self.__CheckValues(entity.GetValue(input_var) * 2, entity.GetValue(output_var))

            # modify the read tensor data
            read_tensor_adaptor.data *= 6

            read_tensor_adaptor.StoreData()
            for entity in container:
                self.__CheckValues(entity.GetValue(input_var), entity.GetValue(output_var) * 3)

    def __CheckValues(self, value_1, value_2):
        if isinstance(value_1, float):
            self.assertEqual(value_1, value_2)
        elif isinstance(value_1, (Kratos.Array3, Kratos.Array4, Kratos.Array6, Kratos.Array9, Kratos.Vector)):
            self.assertVectorAlmostEqual(value_1, value_2)
        elif isinstance(value_1, Kratos.Matrix):
            self.assertMatrixAlmostEqual(value_1, value_2)
        else:
            raise RuntimeError(f"Unsupported types provided [ value_1 = {value_1}, value_2 = {value_2}]")

if __name__ == '__main__':
    KratosUnittest.main()
