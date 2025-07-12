
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
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NODAL_VAUX)

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
            cls.model_part.CreateNewElement("EdgeBasedGradientRecoveryElement3D2N", i + 1, [i + 1, ((i + 1) % 10 + 1)], cls.model_part.GetProperties(i + 1))

        Kratos.VariableUtils().AddDof(Kratos.NODAL_VAUX_X, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.NODAL_VAUX_Y, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.NODAL_VAUX_Z, cls.model_part)

        builder_and_solver = Kratos.ResidualBasedBlockBuilderAndSolver(Kratos.LinearSolverFactory().Create(Kratos.Parameters("""{"solver_type": "skyline_lu_factorization"}""")))
        builder_and_solver.SetUpDofSet(Kratos.ResidualBasedBDFDisplacementScheme(), cls.model_part)
        builder_and_solver.SetUpSystem(cls.model_part)

    def setUp(self):
        Kratos.VariableUtils().ClearNonHistoricalData(self.model_part.Nodes)
        Kratos.VariableUtils().ClearNonHistoricalData(self.model_part.Conditions)
        Kratos.VariableUtils().ClearNonHistoricalData(self.model_part.Elements)
        Kratos.VariableUtils().ClearNonHistoricalData(self.model_part.Properties)

        def value_setter(container, setter_method):
            for entity in container:
                setter_method(entity, Kratos.PRESSURE, entity.Id + 1)
                setter_method(entity, Kratos.VELOCITY, Kratos.Vector([entity.Id * 2 + 1, entity.Id * 3 + 1, entity.Id * 4 + 1]))
                setter_method(entity, Kratos.EXTERNAL_FORCES_VECTOR, Kratos.Vector([entity.Id * 2 + 1, entity.Id * 3 + 1, entity.Id * 4 + 1, entity.Id * 2 + 1, entity.Id * 3 + 1]))
                setter_method(entity, Kratos.NORMAL_SHAPE_DERIVATIVE, Kratos.Matrix([[entity.Id * 2 + 1, entity.Id * 3 + 1], [entity.Id * 4 + 1, entity.Id * 2 + 1], [entity.Id * 3 + 1, entity.Id * 4 + 1]]))

        def flag_setter(container):
            for entity in container:
                entity.Set(Kratos.SLIP, entity.Id % 2)
                entity.Set(Kratos.SELECTED, False)

        value_setter(self.model_part.Nodes, lambda x, y, z: x.SetSolutionStepValue(y, 0, z * 1.5))
        value_setter(self.model_part.Nodes, lambda x, y, z: x.SetSolutionStepValue(y, 1, z * 1.5))
        value_setter(self.model_part.Nodes, lambda x, y, z: x.SetValue(y, z * 2.1))
        value_setter(self.model_part.Properties, lambda x, y, z: x.SetValue(y, z * 6.1))
        value_setter(self.model_part.Conditions, lambda x, y, z: x.SetValue(y, z * 3.1))
        value_setter(self.model_part.Elements, lambda x, y, z: x.SetValue(y, z * 4.1))

        flag_setter(self.model_part.Nodes)
        flag_setter(self.model_part.Conditions)
        flag_setter(self.model_part.Elements)

    def test_ShapeHistoricalVariableTensorAdaptor(self):
        self.__TestTensorAdaptorShape(self.model_part.Nodes, Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor)

    def test_ShapeVariableTensorAdaptorNodes(self):
        self.__TestTensorAdaptorShape(self.model_part.Nodes, Kratos.TensorAdaptors.VariableTensorAdaptor)

    def test_ShapeVariableTensorAdaptorConditions(self):
        self.__TestTensorAdaptorShape(self.model_part.Conditions, Kratos.TensorAdaptors.VariableTensorAdaptor)

    def test_ShapeVariableTensorAdaptorElements(self):
        self.__TestTensorAdaptorShape(self.model_part.Elements, Kratos.TensorAdaptors.VariableTensorAdaptor)

    def test_ShapeVariableTensorAdaptorProperties(self):
        self.__TestTensorAdaptorShape(self.model_part.Properties, Kratos.TensorAdaptors.VariableTensorAdaptor)

    def test_CustomHistoricalVariableTensorAdaptor(self):
        self.__TestCustomShape(self.model_part.Nodes, Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor)

    def test_CustomVariableTensorAdaptorNodes(self):
        self.__TestCustomShape(self.model_part.Nodes, Kratos.TensorAdaptors.VariableTensorAdaptor)

    def test_CustomVariableTensorAdaptorConditions(self):
        self.__TestCustomShape(self.model_part.Conditions, Kratos.TensorAdaptors.VariableTensorAdaptor)

    def test_CustomVariableTensorAdaptorElements(self):
        self.__TestCustomShape(self.model_part.Elements, Kratos.TensorAdaptors.VariableTensorAdaptor)

    def test_CustomVariableTensorAdaptorProperties(self):
        self.__TestCustomShape(self.model_part.Properties, Kratos.TensorAdaptors.VariableTensorAdaptor)

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

        for i, numpy_dtype in enumerate([bool, numpy.uint8, numpy.uint16, numpy.uint32, numpy.uint64, numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64, numpy.float128]):
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

    def test_NumpyArrayAssignment(self):
        t_adaptor_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        t_adaptor_1.CollectData()

        numpy_ta_1 = t_adaptor_1.data

        t_adaptor_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.DENSITY)
        t_adaptor_2.data = numpy_ta_1
        t_adaptor_2.StoreData()

        for node in self.model_part.Nodes:
            self.assertEqual(node.GetValue(Kratos.PRESSURE), node.GetValue(Kratos.DENSITY))

        # now change the values of pressure in every node
        t_adaptor_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        t_adaptor_1.CollectData()
        t_adaptor_1.data *= 2
        t_adaptor_1.StoreData()

        numpy_ta_1 = t_adaptor_1.data
        numpy_ta_2 = t_adaptor_2.data

        # now assign numpy_ta_1 to numpy_ta_2, just only copying the data
        numpy_ta_2[:] = numpy_ta_1

        t_adaptor_2.StoreData()
        for node in self.model_part.Nodes:
            self.assertEqual(node.GetValue(Kratos.PRESSURE), node.GetValue(Kratos.DENSITY))

    def test_TensorAdaptorMoveData(self):
        t_adaptor_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        t_adaptor_1.CollectData()

        numpy_data = t_adaptor_1.MoveData()
        for i, node in enumerate(self.model_part.Nodes):
            self.assertEqual(numpy_data[i], node.GetValue(Kratos.PRESSURE))

        # now the tensor adaptor should be unusable
        with self.assertRaises(RuntimeError):
            t_adaptor_1.data

    def test_NodeVariableTensorAdaptor(self):
        self.__TestVariableTensorAdaptor(self.model_part.Nodes)

    def test_ConditionVariableTensorAdaptor(self):
        self.__TestVariableTensorAdaptor(self.model_part.Conditions)

    def test_ElementVariableTensorAdaptor(self):
        self.__TestVariableTensorAdaptor(self.model_part.Elements)

    def test_PropertyVariableTensorAdaptor(self):
        self.__TestVariableTensorAdaptor(self.model_part.Properties)

    def test_ElementEquationIdsTensorAdaptor(self):
        self.__TestEquationIdsTensorAdaptor(self.model_part.Elements)

    def test_NodalFlagsTensorAdaptor(self):
        self.__TestFlagsTensorAdaptor(self.model_part.Nodes)

    def test_ConditionFlagsTensorAdaptor(self):
        self.__TestFlagsTensorAdaptor(self.model_part.Conditions)

    def test_ElementFlagsTensorAdaptor(self):
        self.__TestFlagsTensorAdaptor(self.model_part.Elements)

    def test_NodePositionTensorAdaptor1(self):
        node_position_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(self.model_part.Nodes, Kratos.Configuration.Current)
        node_position_ta.CollectData()

        numpy_data = node_position_ta.data
        for i, node in enumerate(self.model_part.Nodes):
            self.assertVectorAlmostEqual(numpy_data[i, :], [node.X, node.Y, node.Z])

        node_position_ta_1 = Kratos.TensorAdaptors.NodePositionTensorAdaptor(self.model_part.Nodes, Kratos.Configuration.Current, [2])
        node_position_ta_1.CollectData()
        numpy_data_1 = node_position_ta_1.data
        for i, node in enumerate(self.model_part.Nodes):
            self.assertVectorAlmostEqual(numpy_data_1[i, :], [node.X, node.Y])

        numpy_data_1[:] = numpy_data[:, 0:2] * 2
        node_position_ta_1.StoreData()
        for i, node in enumerate(self.model_part.Nodes):
            self.assertEqual(node.X, numpy_data[i, 0] * 2)
            self.assertEqual(node.Y, numpy_data[i, 1] * 2)
            self.assertEqual(node.Z, numpy_data[i, 2])

    def test_NodePositionTensorAdaptor2(self):
        node_position_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(self.model_part.Nodes, Kratos.Configuration.Initial)
        node_position_ta.CollectData()

        numpy_data = node_position_ta.data
        for i, node in enumerate(self.model_part.Nodes):
            self.assertVectorAlmostEqual(numpy_data[i, :], [node.X0, node.Y0, node.Z0])

        node_position_ta_1 = Kratos.TensorAdaptors.NodePositionTensorAdaptor(self.model_part.Nodes, Kratos.Configuration.Initial, [2])
        node_position_ta_1.CollectData()
        numpy_data_1 = node_position_ta_1.data
        for i, node in enumerate(self.model_part.Nodes):
            self.assertVectorAlmostEqual(numpy_data_1[i, :], [node.X0, node.Y0])

        numpy_data_1[:] = numpy_data[:, 0:2] * 2
        node_position_ta_1.StoreData()
        for i, node in enumerate(self.model_part.Nodes):
            self.assertEqual(node.X0, numpy_data[i, 0] * 2)
            self.assertEqual(node.Y0, numpy_data[i, 1] * 2)
            self.assertEqual(node.Z0, numpy_data[i, 2])

    def __TestFlagsTensorAdaptor(self, container):
        tensor_adaptor_read = Kratos.TensorAdaptors.FlagsTensorAdaptor(container, Kratos.SLIP)
        tensor_adaptor_read.CollectData()

        numpy_data_read = tensor_adaptor_read.data
        for i, entity in enumerate(container):
            self.assertEqual(numpy_data_read[i], entity.Is(Kratos.SLIP))

        tensor_adaptor_write = Kratos.TensorAdaptors.FlagsTensorAdaptor(tensor_adaptor_read.GetContainer(), Kratos.SELECTED)
        numpy_data_write = tensor_adaptor_write.data

        numpy_data_write[:] = numpy.invert(numpy_data_read)
        tensor_adaptor_write.StoreData()

        for entity in container:
            self.assertEqual(entity.Is(Kratos.SLIP), not entity.Is(Kratos.SELECTED))

    def __TestEquationIdsTensorAdaptor(self, container):
        tensor_adaptor = Kratos.TensorAdaptors.EquationIdsTensorAdaptor(container, self.model_part.ProcessInfo)
        tensor_adaptor.CollectData()

        numpy_data = tensor_adaptor.data

        for i, entity in enumerate(container):
            self.assertVectorAlmostEqual(numpy_data[i, :], entity.EquationIdVector(self.model_part.ProcessInfo))

        with self.assertRaises(RuntimeError):
            tensor_adaptor.StoreData()

    def __TestVariableTensorAdaptor(self, container):
        for input_var, output_var in zip(self.input_list_of_variables, self.output_list_of_variables):
            read_tensor_adaptor = Kratos.TensorAdaptors.VariableTensorAdaptor(container, input_var)
            read_tensor_adaptor.CollectData()

            write_tensor_adaptor = Kratos.TensorAdaptors.VariableTensorAdaptor(read_tensor_adaptor.GetContainer(), output_var, data_shape=read_tensor_adaptor.DataShape())

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

    def __TestTensorAdaptorShape(self, container, TensorAdaptorType):
        t_adaptor_1 = TensorAdaptorType(container, Kratos.VELOCITY)
        self.assertVectorAlmostEqual(t_adaptor_1.Shape(), [len(container), 3])
        self.assertVectorAlmostEqual(t_adaptor_1.DataShape(), [3])

        t_adaptor_1 = TensorAdaptorType(container, Kratos.NORMAL_SHAPE_DERIVATIVE)
        self.assertVectorAlmostEqual(t_adaptor_1.Shape(), [len(container), 3, 2])
        self.assertVectorAlmostEqual(t_adaptor_1.DataShape(), [3, 2])

        t_adaptor_1 = TensorAdaptorType(container, Kratos.PK2_STRESS_TENSOR)
        self.assertVectorAlmostEqual(t_adaptor_1.Shape(), [len(container), 0, 0])
        self.assertVectorAlmostEqual(t_adaptor_1.DataShape(), [0, 0])

        t_adaptor_1 = TensorAdaptorType(container, Kratos.PK2_STRESS_TENSOR, data_shape=[5,6])
        self.assertVectorAlmostEqual(t_adaptor_1.Shape(), [len(container), 5, 6])
        self.assertVectorAlmostEqual(t_adaptor_1.DataShape(), [5, 6])

    def __TestCustomShape(self, container, TensorAdaptorType):
        t_adaptor_1 = TensorAdaptorType(container, Kratos.VELOCITY)
        t_adaptor_1.CollectData()

        t_adaptor_2 = TensorAdaptorType(container, Kratos.VELOCITY, data_shape=[2])
        t_adaptor_2.CollectData()

        t_adaptor_3 = TensorAdaptorType(container, Kratos.VELOCITY, data_shape=[1])
        t_adaptor_3.CollectData()

        t_adaptor_4 = TensorAdaptorType(container, Kratos.VELOCITY_X)
        t_adaptor_4.CollectData()

        self.assertVectorAlmostEqual(t_adaptor_4.data, t_adaptor_3.data)
        self.assertVectorAlmostEqual(t_adaptor_3.data, t_adaptor_1.data[:, 0])
        for v1, v2 in zip(t_adaptor_2.data, t_adaptor_1.data[:, :2]):
            self.assertVectorAlmostEqual(v1, v2)

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
