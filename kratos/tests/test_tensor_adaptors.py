
import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as KratosUtilities
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestTensorAdaptors(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.SetBufferSize(2)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NODAL_VAUX)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)

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
                setter_method(entity, Kratos.NORMAL_SHAPE_DERIVATIVE, Kratos.Matrix([[entity.Id * 2 + 1, entity.Id * 3 + 1, entity.Id, 2 * entity.Id], [entity.Id * 4 + 1, entity.Id * 2 + 1, entity.Id, 3 * entity.Id], [entity.Id * 3 + 1, entity.Id * 4 + 1, entity.Id, 4 * entity.Id]]))

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

        with self.assertRaises(RuntimeError):
            t_adaptor_1.Check()

        t_adaptor_2.StoreData()

        with self.assertRaises(RuntimeError):
            t_adaptor_3 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.CONSTITUTIVE_MATRIX)

        t_adaptor_3 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.CONSTITUTIVE_MATRIX, data_shape=[2,3])

        with self.assertRaises(RuntimeError):
            t_adaptor_3.Check()

        with self.assertRaises(RuntimeError):
            t_adaptor_5 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.NORMAL_SHAPE_DERIVATIVE, step_index=6)

        t_adaptor_5 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.NORMAL_SHAPE_DERIVATIVE, data_shape=[2,3], step_index=6)
        with self.assertRaises(RuntimeError):
            t_adaptor_5.Check()

        with self.assertRaises(RuntimeError):
            Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Conditions, Kratos.VELOCITY), Kratos.PRESSURE)

        with self.assertRaises(RuntimeError):
            Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.VELOCITY), Kratos.PRESSURE)

        with self.assertRaises(RuntimeError):
            Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Properties, Kratos.VELOCITY), Kratos.PRESSURE)

    def test_SupportedNdArrays(self):
        t_adaptor = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.DENSITY)

        for i, numpy_dtype in enumerate([bool, numpy.uint8, numpy.uint16, numpy.uint32, numpy.uint64, numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64]):
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

    def test_NodalVariableTensorAdaptorPartial(self):
        self.__TestVariableTensorAdaptorPartial(self.model_part.Nodes)

    def test_ConditionVariableTensorAdaptorPartial(self):
        self.__TestVariableTensorAdaptorPartial(self.model_part.Conditions)

    def test_ElementVariableTensorAdaptorPartial(self):
        self.__TestVariableTensorAdaptorPartial(self.model_part.Elements)

    def test_PropertiesVariableTensorAdaptorPartial(self):
        self.__TestVariableTensorAdaptorPartial(self.model_part.Properties)

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

        with self.assertRaises(RuntimeError):
            Kratos.TensorAdaptors.NodePositionTensorAdaptor(Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Conditions, Kratos.VELOCITY), Kratos.Configuration.Initial)

        with self.assertRaises(RuntimeError):
            Kratos.TensorAdaptors.NodePositionTensorAdaptor(Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.VELOCITY), Kratos.Configuration.Initial)

        with self.assertRaises(RuntimeError):
            Kratos.TensorAdaptors.NodePositionTensorAdaptor(Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Properties, Kratos.VELOCITY), Kratos.Configuration.Initial)

    def test_BaseTensorAdaptor(self):
        var_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)
        var_ta.Check()
        var_ta.CollectData()

        copied_ta = Kratos.TensorAdaptors.DoubleTensorAdaptor(var_ta)
        base_ta = Kratos.TensorAdaptors.DoubleTensorAdaptor(var_ta, copy=False)
        base_ta.data += 1.0

        var_ta.StoreData()

        for i, node in enumerate(self.model_part.Nodes):
            self.assertEqual(node.GetValue(Kratos.VELOCITY_X), copied_ta.data[i, 0] + 1)
            self.assertEqual(node.GetValue(Kratos.VELOCITY_Y), copied_ta.data[i, 1] + 1)
            self.assertEqual(node.GetValue(Kratos.VELOCITY_Z), copied_ta.data[i, 2] + 1)

    def test_CopyVariableTensorAdaptor(self):
        self.__TestCopyTensorAdaptor(Kratos.TensorAdaptors.VariableTensorAdaptor, lambda x, y: x.GetValue(y))

    def test_CopyHistoricalVariableTensorAdaptor(self):
        self.__TestCopyTensorAdaptor(Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor, lambda x, y: x.GetSolutionStepValue(y))

    def test_FStyleNodalPositionTensorAdaptor(self):
        ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(self.model_part.Nodes, Kratos.Configuration.Initial)
        ta.CollectData()

        numpy_data = ta.data
        fortran_numpy_data = numpy.asfortranarray(numpy_data)
        with self.assertRaises(RuntimeError):
            ta.data = fortran_numpy_data

    def test_TensorAdaptorWithNDData(self):
        np_array = numpy.ones((self.model_part.NumberOfNodes(), 3))
        dd_array = Kratos.DoubleNDData(np_array, copy=False)
        ta_no_copy = Kratos.TensorAdaptors.DoubleTensorAdaptor(self.model_part.Nodes, dd_array, copy = False)
        ta_copy = Kratos.TensorAdaptors.DoubleTensorAdaptor(self.model_part.Nodes, dd_array)

        self.assertEqual(numpy.linalg.norm(ta_no_copy.data - np_array), 0.0)
        self.assertEqual(numpy.linalg.norm(ta_copy.data - np_array), 0.0)
        np_array[0, 0] = -500.0
        self.assertEqual(numpy.linalg.norm(ta_no_copy.data - np_array), 0.0)
        self.assertEqual(numpy.sum(ta_copy.data - np_array), 501)

        Kratos.TensorAdaptors.NodePositionTensorAdaptor(ta_copy, Kratos.Configuration.Initial, copy=False).StoreData()
        for node in self.model_part.Nodes:
            self.assertEqual(node.X0, 1.0)
            self.assertEqual(node.Y0, 1.0)
            self.assertEqual(node.Z0, 1.0)

    def test_CombinedTensorAdaptorAxis1(self):
        x_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(self.model_part.Nodes, Kratos.Configuration.Initial)
        p_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        u_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)

        combined_ta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([p_ta, x_ta, u_ta], axis=1)
        combined_ta.Check()
        combined_ta.CollectData()

        numpy_concat = numpy.concatenate([p_ta.data.reshape((len(self.model_part.Nodes), 1)), x_ta.data, u_ta.data], axis=1)

        self.assertVectorAlmostEqual(numpy_concat.shape, combined_ta.Shape())

        for i in range(self.model_part.NumberOfNodes()):
            for j in range(7):
                self.assertEqual(numpy_concat[i, j], combined_ta.data[i, j])

        self.assertFalse(combined_ta.HasContainer())
        with self.assertRaises(RuntimeError):
            combined_ta.GetContainer()

        combined_ta.data *= 1.5
        combined_ta.StoreData()

        self.assertEqual(numpy.linalg.norm(p_ta.data - numpy_concat[:, 0  ] * 1.5), 0.0)
        self.assertEqual(numpy.linalg.norm(x_ta.data - numpy_concat[:, 1:4] * 1.5), 0.0)
        self.assertEqual(numpy.linalg.norm(u_ta.data - numpy_concat[:, 4:7] * 1.5), 0.0)

        for i, node in enumerate(self.model_part.Nodes):
            self.assertEqual(node.GetValue(Kratos.PRESSURE), numpy_concat[i, 0] * 1.5)
            self.assertEqual(node.X0, numpy_concat[i, 1] * 1.5)
            self.assertEqual(node.Y0, numpy_concat[i, 2] * 1.5)
            self.assertEqual(node.Z0, numpy_concat[i, 3] * 1.5)
            self.assertEqual(node.GetValue(Kratos.VELOCITY_X), numpy_concat[i, 4] * 1.5)
            self.assertEqual(node.GetValue(Kratos.VELOCITY_Y), numpy_concat[i, 5] * 1.5)
            self.assertEqual(node.GetValue(Kratos.VELOCITY_Z), numpy_concat[i, 6] * 1.5)

        # check for recursive collect and store data
        x_ta.data *= 2.0
        u_ta.data *= 2.0
        p_ta.data *= 2.0
        x_ta.StoreData()
        u_ta.StoreData()
        p_ta.StoreData()
        x_ta.data /= 2.0
        u_ta.data /= 2.0
        p_ta.data /= 2.0

        copy_combined_ta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(combined_ta, perform_collect_data_recursively=False, perform_store_data_recursively=False)
        copy_combined_ta.CollectData()
        for i, node in enumerate(self.model_part.Nodes):
            self.assertAlmostEqual(node.GetValue(Kratos.PRESSURE), copy_combined_ta.data[i, 0] * 2.0)
            self.assertAlmostEqual(node.X0, copy_combined_ta.data[i, 1] * 2.0)
            self.assertAlmostEqual(node.Y0, copy_combined_ta.data[i, 2] * 2.0)
            self.assertAlmostEqual(node.Z0, copy_combined_ta.data[i, 3] * 2.0)
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_X), copy_combined_ta.data[i, 4] * 2.0)
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_Y), copy_combined_ta.data[i, 5] * 2.0)
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_Z), copy_combined_ta.data[i, 6] * 2.0)

        copy_combined_ta.data *= 3.0
        copy_combined_ta.StoreData()
        for i, node in enumerate(self.model_part.Nodes):
            self.assertAlmostEqual(node.GetValue(Kratos.PRESSURE), p_ta.data[i] * 2.0 / 3.0)
            self.assertAlmostEqual(node.X0, x_ta.data[i, 0] * 2.0 / 3.0)
            self.assertAlmostEqual(node.Y0, x_ta.data[i, 1] * 2.0 / 3.0)
            self.assertAlmostEqual(node.Z0, x_ta.data[i, 2] * 2.0 / 3.0)
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_X), u_ta.data[i, 0] * 2.0 / 3.0)
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_Y), u_ta.data[i, 1] * 2.0 / 3.0)
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_Z), u_ta.data[i, 2] * 2.0 / 3.0)

    def test_CombinedTensorAdaptorAxis0(self):
        u_nodes_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)
        u_elems_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.VELOCITY)

        combined_ta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([u_nodes_ta, u_elems_ta], axis=0)
        combined_ta.Check()
        combined_ta.CollectData()

        numpy_concat = numpy.concatenate([u_nodes_ta.data, u_elems_ta.data])

        self.assertVectorAlmostEqual(numpy_concat.shape, combined_ta.Shape())

        for i in range(self.model_part.NumberOfNodes() + self.model_part.NumberOfElements()):
            for j in range(3):
                self.assertEqual(numpy_concat[i, j], combined_ta.data[i, j])

        self.assertFalse(combined_ta.HasContainer())
        with self.assertRaises(RuntimeError):
            combined_ta.GetContainer()

        combined_ta.data *= 2.3
        combined_ta.StoreData()
        for i, node in enumerate(self.model_part.Nodes):
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_X), numpy_concat[i, 0] * 2.3)
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_Y), numpy_concat[i, 1] * 2.3)
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_Z), numpy_concat[i, 2] * 2.3)

        for i, node in enumerate(self.model_part.Elements):
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_X), numpy_concat[i + 10, 0] * 2.3)
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_Y), numpy_concat[i + 10, 1] * 2.3)
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_Z), numpy_concat[i + 10, 2] * 2.3)

    def test_CombinedTensorAdaptorRavel(self):
        u_nodes_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)
        u_elems_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.PRESSURE)

        combined_ta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([u_nodes_ta, u_elems_ta])
        combined_ta.Check()
        combined_ta.CollectData()

        numpy_concat = numpy.concatenate([u_nodes_ta.data.ravel(), u_elems_ta.data.ravel()])

        self.assertVectorAlmostEqual(numpy_concat.shape, combined_ta.Shape())

        for i in range(numpy_concat.shape[0]):
            self.assertEqual(numpy_concat[i], combined_ta.data[i])

        self.assertFalse(combined_ta.HasContainer())
        with self.assertRaises(RuntimeError):
            combined_ta.GetContainer()

        combined_ta.data *= 2.8
        combined_ta.StoreData()
        index = 0
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_X), numpy_concat[index] * 2.8)
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_Y), numpy_concat[index + 1] * 2.8)
            self.assertAlmostEqual(node.GetValue(Kratos.VELOCITY_Z), numpy_concat[index + 2] * 2.8)
            index += 3

        for i, element in enumerate(self.model_part.Elements):
            self.assertAlmostEqual(element.GetValue(Kratos.PRESSURE), numpy_concat[30 + i] * 2.8)

    def test_CombinedTensorAdaptorNested(self):
        p_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        hist_p_ta = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        u_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)
        hist_u_ta = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)

        combined_ta_1 = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([hist_p_ta, u_ta], axis=1)
        combined_ta_2 = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([p_ta, hist_u_ta], axis=1)
        combined_ta_3 = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([combined_ta_1, combined_ta_2], axis=1)

        combined_ta_3.Check()
        combined_ta_3.CollectData()

        numpy_data = numpy.array(combined_ta_3.data)

        for i, node in enumerate(self.model_part.Nodes):
            self.assertEqual(node.GetSolutionStepValue(Kratos.PRESSURE), combined_ta_3.data[i, 0])
            self.assertEqual(node.GetValue(Kratos.VELOCITY_X), combined_ta_3.data[i, 1])
            self.assertEqual(node.GetValue(Kratos.VELOCITY_Y), combined_ta_3.data[i, 2])
            self.assertEqual(node.GetValue(Kratos.VELOCITY_Z), combined_ta_3.data[i, 3])
            self.assertEqual(node.GetValue(Kratos.PRESSURE), combined_ta_3.data[i, 4])
            self.assertEqual(node.GetSolutionStepValue(Kratos.VELOCITY_X), combined_ta_3.data[i, 5])
            self.assertEqual(node.GetSolutionStepValue(Kratos.VELOCITY_Y), combined_ta_3.data[i, 6])
            self.assertEqual(node.GetSolutionStepValue(Kratos.VELOCITY_Z), combined_ta_3.data[i, 7])

        combined_ta_3.data *= 4.5
        combined_ta_3.StoreData()

        for i, node in enumerate(self.model_part.Nodes):
            self.assertEqual(node.GetSolutionStepValue(Kratos.PRESSURE), numpy_data[i, 0] * 4.5)
            self.assertEqual(node.GetValue(Kratos.VELOCITY_X), numpy_data[i, 1] * 4.5)
            self.assertEqual(node.GetValue(Kratos.VELOCITY_Y), numpy_data[i, 2] * 4.5)
            self.assertEqual(node.GetValue(Kratos.VELOCITY_Z), numpy_data[i, 3] * 4.5)
            self.assertEqual(node.GetValue(Kratos.PRESSURE), numpy_data[i, 4] * 4.5)
            self.assertEqual(node.GetSolutionStepValue(Kratos.VELOCITY_X), numpy_data[i, 5] * 4.5)
            self.assertEqual(node.GetSolutionStepValue(Kratos.VELOCITY_Y), numpy_data[i, 6] * 4.5)
            self.assertEqual(node.GetSolutionStepValue(Kratos.VELOCITY_Z), numpy_data[i, 7] * 4.5)

    def test_CombinedTensorAdaptorGetTensorAdaptors(self):
        u_nodes_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)
        u_elems_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.PRESSURE)

        combined_ta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([u_nodes_ta, u_elems_ta])

        self.assertEqual(id(u_nodes_ta), id(combined_ta.GetTensorAdaptors()[0]))
        self.assertEqual(id(u_elems_ta), id(combined_ta.GetTensorAdaptors()[1]))

    def test_CombinedTensorAdaptorCopy(self):
        u_nodes_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY)
        u_elems_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.PRESSURE)

        combined_ta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([u_nodes_ta, u_elems_ta])
        combined_ta.Check()
        combined_ta.CollectData()

        copied_combined_ta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(combined_ta)
        copied_combined_ta.CollectData()

        copied_combined_ta.data *= 2.0

        self.assertAlmostEqual(numpy.linalg.norm(combined_ta.data * 2.0 - copied_combined_ta.data), 0.0)

        ref_combined_ta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(combined_ta, copy=False)
        ref_combined_ta.data *= 3.0
        self.assertAlmostEqual(numpy.linalg.norm(ref_combined_ta.data - combined_ta.data), 0.0)

    def test_NodalNeighboursTensorAdaptorCondition(self):
        ta = Kratos.TensorAdaptors.NodalNeighbourCountTensorAdaptor(self.model_part.Nodes, self.model_part.Conditions)
        ta.Check()
        ta.CollectData()

        for i, v in enumerate(ta.data):
            if (i != 0 or i != ta.data.size - 1):
                self.assertEqual(v, 2)
            else:
                self.assertEqual(v, 1)

        with self.assertRaises(RuntimeError):
            ta.StoreData()

    def test_NodalNeighboursTensorAdaptorElement(self):
        ta = Kratos.TensorAdaptors.NodalNeighbourCountTensorAdaptor(self.model_part.Nodes, self.model_part.Elements)
        ta.Check()
        ta.CollectData()

        for i, v in enumerate(ta.data):
            if (i != 0 or i != ta.data.size - 1):
                self.assertEqual(v, 2)
            else:
                self.assertEqual(v, 1)

        with self.assertRaises(RuntimeError):
            ta.StoreData()

    def testGeometryMetricsTensorAdaptorDomainSizeCondition(self):
        ta = Kratos.TensorAdaptors.GeometryMetricsTensorAdaptor(self.model_part.Conditions, Kratos.TensorAdaptors.GeometryMetricsTensorAdaptor.DomainSize)
        ta.CollectData()
        for i, condition in enumerate(self.model_part.Conditions):
            self.assertEqual(ta.data[i], condition.GetGeometry().DomainSize())

    def testGeometryMetricsTensorAdaptorDomainSizeElement(self):
        ta = Kratos.TensorAdaptors.GeometryMetricsTensorAdaptor(self.model_part.Elements, Kratos.TensorAdaptors.GeometryMetricsTensorAdaptor.DomainSize)
        ta.CollectData()
        for i, element in enumerate(self.model_part.Elements):
            self.assertEqual(ta.data[i], element.GetGeometry().DomainSize())

    def testGeometryMetricsTensorAdaptorDomainSizeCondition_Empty(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        ta = Kratos.TensorAdaptors.GeometryMetricsTensorAdaptor(model_part.Conditions, Kratos.TensorAdaptors.GeometryMetricsTensorAdaptor.DomainSize)
        ta.CollectData()
        self.assertEqual(ta.data.shape, (0,))

    def testGeometryMetricsTensorAdaptorDomainSizeElement_Empty(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        ta = Kratos.TensorAdaptors.GeometryMetricsTensorAdaptor(model_part.Elements, Kratos.TensorAdaptors.GeometryMetricsTensorAdaptor.DomainSize)
        ta.CollectData()
        self.assertEqual(ta.data.shape, (0,))

    def __TestCopyTensorAdaptor(self, tensor_adaptor_type, value_getter):
        var_ta_orig = tensor_adaptor_type(self.model_part.Nodes, Kratos.VELOCITY, data_shape=[2])
        var_ta_orig.Check()
        var_ta_orig.CollectData()

        var_ta_dest = tensor_adaptor_type(var_ta_orig, Kratos.ACCELERATION, copy=False)
        var_ta_dest.data += 1.0
        var_ta_dest.StoreData()

        for node in self.model_part.Nodes:
            self.assertEqual(value_getter(node, Kratos.ACCELERATION_X), value_getter(node, Kratos.VELOCITY_X) + 1.0)
            self.assertEqual(value_getter(node, Kratos.ACCELERATION_Y), value_getter(node, Kratos.VELOCITY_Y) + 1.0)
            self.assertEqual(value_getter(node, Kratos.ACCELERATION_Z), 0.0)

        # now calling the dest CollectData of var_ta_dest. Since we did not state copy = true,
        # it should change the storage of the var_ta_orig as well.
        var_ta_dest.Check()
        var_ta_dest.CollectData()
        var_ta_orig.StoreData()
        for node in self.model_part.Nodes:
            self.assertEqual(value_getter(node, Kratos.ACCELERATION_X), value_getter(node, Kratos.VELOCITY_X))
            self.assertEqual(value_getter(node, Kratos.ACCELERATION_Y), value_getter(node, Kratos.VELOCITY_Y))

        # now checking the copying
        var_ta_copy = tensor_adaptor_type(var_ta_orig, Kratos.DISPLACEMENT) # default behavior is copying
        var_ta_copy.data += 2.0
        var_ta_copy.StoreData()
        for node in self.model_part.Nodes:
            self.assertEqual(value_getter(node, Kratos.DISPLACEMENT_X), value_getter(node, Kratos.VELOCITY_X) + 2.0)
            self.assertEqual(value_getter(node, Kratos.DISPLACEMENT_Y), value_getter(node, Kratos.VELOCITY_Y) + 2.0)
            self.assertEqual(value_getter(node, Kratos.DISPLACEMENT_Z), 0.0)

        # now reading data to the copied. It should not change the var_ta_orig
        var_ta_copy.Check()
        var_ta_copy.CollectData()
        var_ta_orig.StoreData()
        for node in self.model_part.Nodes:
            self.assertEqual(value_getter(node, Kratos.ACCELERATION_X), value_getter(node, Kratos.VELOCITY_X))
            self.assertEqual(value_getter(node, Kratos.ACCELERATION_Y), value_getter(node, Kratos.VELOCITY_Y))

    def __TestFlagsTensorAdaptor(self, container):
        tensor_adaptor_read = Kratos.TensorAdaptors.FlagsTensorAdaptor(container, Kratos.SLIP)
        tensor_adaptor_read.CollectData()

        numpy_data_read = tensor_adaptor_read.data
        for i, entity in enumerate(container):
            self.assertEqual(numpy_data_read[i], entity.Is(Kratos.SLIP))

        tensor_adaptor_write = Kratos.TensorAdaptors.FlagsTensorAdaptor(tensor_adaptor_read.GetContainer(), Kratos.SELECTED)
        numpy_data_write = tensor_adaptor_write.data

        numpy_data_write[:] = numpy.invert(numpy.asarray(numpy_data_read, dtype=bool))
        tensor_adaptor_write.StoreData()

        for entity in container:
            self.assertEqual(entity.Is(Kratos.SLIP), not entity.Is(Kratos.SELECTED))

        # testing undefined status
        for entity in container:
            if entity.Id % 3 == 0:
                entity.Set(Kratos.VISITED, entity.Id % 2)

        tensor_adaptor_undef_read = Kratos.TensorAdaptors.FlagsTensorAdaptor(container, Kratos.VISITED)
        tensor_adaptor_undef_read.CollectData()

        for entity, value in zip(container, tensor_adaptor_undef_read.data):
            self.assertEqual(entity.IsDefined(Kratos.VISITED), value != -1)
            if value > 0:
                self.assertEqual(entity.Is(Kratos.VISITED), bool(value))

        tensor_adaptor_undef_write = Kratos.TensorAdaptors.FlagsTensorAdaptor(tensor_adaptor_undef_read, Kratos.STRUCTURE)
        tensor_adaptor_undef_write.StoreData()
        for entity in container:
            self.assertEqual(entity.IsDefined(Kratos.VISITED), entity.IsDefined(Kratos.STRUCTURE))
            self.assertEqual(entity.Is(Kratos.VISITED), entity.Is(Kratos.STRUCTURE))

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

    def __TestVariableTensorAdaptorPartial(self, container):
        read_ta_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.VELOCITY)
        write_ta_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.ACCELERATION)
        read_ta_1.CollectData()
        write_ta_1.data = read_ta_1.data
        write_ta_1.StoreData()

        read_ta_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.VELOCITY, data_shape=[2])
        read_ta_2.CollectData()
        write_ta_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.ACCELERATION, data_shape=read_ta_2.DataShape())
        write_ta_2.data = read_ta_2.data * 2 # will only modify the first 2 components of the ACCELERATION
        write_ta_2.StoreData()

        for entity in container:
            input = entity.GetValue(Kratos.VELOCITY)
            output = entity.GetValue(Kratos.ACCELERATION)
            self.assertEqual(input[0] * 2, output[0])
            self.assertEqual(input[1] * 2, output[1])
            self.assertEqual(input[2], output[2])

        ## Vector types
        read_ta_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.EXTERNAL_FORCES_VECTOR)
        write_ta_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.INTERNAL_FORCES_VECTOR, data_shape=read_ta_1.DataShape())
        read_ta_1.CollectData()
        write_ta_1.data = read_ta_1.data
        write_ta_1.StoreData()

        read_ta_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.EXTERNAL_FORCES_VECTOR, data_shape=[3])
        read_ta_2.CollectData()
        write_ta_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.INTERNAL_FORCES_VECTOR, data_shape=read_ta_2.DataShape())
        write_ta_2.data = read_ta_2.data * 2 # will only modify the first 3 components of the INTERNAL_FORCES_VECTOR
        write_ta_2.StoreData()

        for entity in container:
            input = entity.GetValue(Kratos.EXTERNAL_FORCES_VECTOR)
            output = entity.GetValue(Kratos.INTERNAL_FORCES_VECTOR)
            self.assertEqual(input[0] * 2, output[0])
            self.assertEqual(input[1] * 2, output[1])
            self.assertEqual(input[2] * 2, output[2])
            self.assertEqual(input[3], output[3])
            self.assertEqual(input[4], output[4])

        write_ta_3 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.BDF_COEFFICIENTS, data_shape=read_ta_2.DataShape())
        write_ta_3.data = read_ta_2.data * 2 # will only have the first 3 components of the BDF_COEFFICIENTS
        write_ta_3.StoreData()

        for entity in container:
            input = entity.GetValue(Kratos.EXTERNAL_FORCES_VECTOR)
            output = entity.GetValue(Kratos.BDF_COEFFICIENTS)
            self.assertEqual(output.Size(), 3)
            self.assertEqual(input.Size(), 5)
            self.assertEqual(input[0] * 2, output[0])
            self.assertEqual(input[1] * 2, output[1])
            self.assertEqual(input[2] * 2, output[2])

        # Matrix type
        read_ta_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.NORMAL_SHAPE_DERIVATIVE)
        write_ta_1 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.PK2_STRESS_TENSOR, data_shape=read_ta_1.DataShape())
        read_ta_1.CollectData()
        write_ta_1.data = read_ta_1.data
        write_ta_1.StoreData()

        read_ta_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.NORMAL_SHAPE_DERIVATIVE, data_shape=[2, 3])
        read_ta_2.CollectData()
        write_ta_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.PK2_STRESS_TENSOR, data_shape=read_ta_2.DataShape())
        write_ta_2.data = read_ta_2.data * 3 # will only modify the first (2,3) components of the PK2_STRESS_TENSOR
        write_ta_2.StoreData()

        for entity in container:
            input = entity.GetValue(Kratos.NORMAL_SHAPE_DERIVATIVE)
            output = entity.GetValue(Kratos.PK2_STRESS_TENSOR)
            self.assertEqual(input.Size1(), 3)
            self.assertEqual(input.Size2(), 4)
            self.assertEqual(output.Size1(), 3)
            self.assertEqual(output.Size2(), 4)
            self.assertEqual(input[0, 0] * 3, output[0, 0])
            self.assertEqual(input[0, 1] * 3, output[0, 1])
            self.assertEqual(input[0, 2] * 3, output[0, 2])
            self.assertEqual(input[0, 3], output[0, 3])
            self.assertEqual(input[1, 0] * 3, output[1, 0])
            self.assertEqual(input[1, 1] * 3, output[1, 1])
            self.assertEqual(input[1, 2] * 3, output[1, 2])
            self.assertEqual(input[1, 3], output[1, 3])
            self.assertEqual(input[2, 0], output[2, 0])
            self.assertEqual(input[2, 1], output[2, 1])
            self.assertEqual(input[2, 2], output[2, 2])
            self.assertEqual(input[2, 3], output[2, 3])

        write_ta_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(container, Kratos.CAUCHY_STRESS_TENSOR, data_shape=read_ta_2.DataShape())
        write_ta_2.data = read_ta_2.data * 3 # will only have the first (2,3) components of the PK2_STRESS_TENSOR
        write_ta_2.StoreData()

        for entity in container:
            input = entity.GetValue(Kratos.NORMAL_SHAPE_DERIVATIVE)
            output = entity.GetValue(Kratos.CAUCHY_STRESS_TENSOR)
            self.assertEqual(input.Size1(), 3)
            self.assertEqual(input.Size2(), 4)
            self.assertEqual(output.Size1(), 2)
            self.assertEqual(output.Size2(), 3)
            self.assertEqual(input[0, 0] * 3, output[0, 0])
            self.assertEqual(input[0, 1] * 3, output[0, 1])
            self.assertEqual(input[0, 2] * 3, output[0, 2])
            self.assertEqual(input[1, 0] * 3, output[1, 0])
            self.assertEqual(input[1, 1] * 3, output[1, 1])
            self.assertEqual(input[1, 2] * 3, output[1, 2])

    def __TestTensorAdaptorShape(self, container, TensorAdaptorType):
        t_adaptor_1 = TensorAdaptorType(container, Kratos.VELOCITY)
        self.assertVectorAlmostEqual(t_adaptor_1.Shape(), [len(container), 3])
        self.assertVectorAlmostEqual(t_adaptor_1.DataShape(), [3])

        t_adaptor_1 = TensorAdaptorType(container, Kratos.NORMAL_SHAPE_DERIVATIVE)
        self.assertVectorAlmostEqual(t_adaptor_1.Shape(), [len(container), 3, 4])
        self.assertVectorAlmostEqual(t_adaptor_1.DataShape(), [3, 4])

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
