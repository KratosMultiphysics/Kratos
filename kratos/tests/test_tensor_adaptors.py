
import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as KratosUtilities
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestTensorAdaptors(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.SetBufferSize(2)

        cls.input_list_of_variables = [
                Kratos.PRESSURE,                # double
                Kratos.VELOCITY,                # array3
                Kratos.EXTERNAL_FORCES_VECTOR,  # vector
                Kratos.NORMAL_SHAPE_DERIVATIVE, # matrix,
                Kratos.INITIAL_STRAIN           # non-initialized variable
            ]

        cls.output_list_of_variables = [
                Kratos.DENSITY,                 # double
                Kratos.ACCELERATION,            # array3
                Kratos.INTERNAL_FORCES_VECTOR,  # vector
                Kratos.PK2_STRESS_TENSOR,       # matrix
                Kratos.RECOVERED_STRESS         # non-initialized variable
            ]

        for var in cls.input_list_of_variables:
            cls.model_part.AddNodalSolutionStepVariable(var)

        for var in cls.output_list_of_variables:
            cls.model_part.AddNodalSolutionStepVariable(var)

        for i in range(10):
            cls.model_part.CreateNewNode(i + 1, i, i + 1, i + 2)

    def setUp(self):
        def setter(container, setter_method):
            for entity in container:
                setter_method(entity, Kratos.PRESSURE, entity.Id + 1)
                setter_method(entity, Kratos.VELOCITY, Kratos.Vector([entity.Id * 2 + 1, entity.Id * 3 + 1, entity.Id * 4 + 1]))
                setter_method(entity, Kratos.EXTERNAL_FORCES_VECTOR, Kratos.Vector([entity.Id * 2 + 1, entity.Id * 3 + 1, entity.Id * 4 + 1, entity.Id * 2 + 1, entity.Id * 3 + 1]))
                setter_method(entity, Kratos.NORMAL_SHAPE_DERIVATIVE, Kratos.Matrix([[entity.Id * 2 + 1, entity.Id * 3 + 1], [entity.Id * 4 + 1, entity.Id * 2 + 1], [entity.Id * 3 + 1, entity.Id * 4 + 1]]))

        setter(self.model_part.Nodes, lambda x, y, z: x.SetSolutionStepValue(y, 0, z))
        setter(self.model_part.Nodes, lambda x, y, z: x.SetSolutionStepValue(y, 1, z))

    def test_NodeHistoricalVariableTensorAdaptorShape(self):
        t = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY, step_index=0, shape=[-1, 3])
        # with self.assertRaises(RuntimeError):
            # Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY, [3, 3], 0)
            # Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.VELOCITY, [3, 3, 3], 0)

    #     self.assertEqual(Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.EXTERNAL_FORCES_VECTOR, 0).Shape(), [10, 5])
    #     self.assertEqual(Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.INTERNAL_FORCES_VECTOR, 0).Shape(), [10, 0])
    #     self.assertEqual(Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.INTERNAL_FORCES_VECTOR, [10, 6], 0).Shape(), [10, 6])

    def test_NodalHistoricalVariableTensorAdaptor1(self):
        self.__TestVariableTensorAdaptor(self.model_part.Nodes, lambda x, y: x.GetSolutionStepValue(y), step_index = 0)

    def test_NodalHistoricalVariableTensorAdaptor2(self):
        self.__TestVariableTensorAdaptor(self.model_part.Nodes, lambda x, y: x.GetSolutionStepValue(y, 1), step_index = 1)

    def __TestVariableTensorAdaptor(self, container, getter_method, **kwargs):
        for input_var, output_var in zip(self.input_list_of_variables, self.output_list_of_variables):
            read_tensor_adaptor = Kratos.TensorAdaptors.VariableTensorAdaptor(container, input_var, **kwargs)

            read_tensor_adaptor.CollectData()

            write_tensor_adaptor = Kratos.TensorAdaptors.VariableTensorAdaptor(read_tensor_adaptor.GetContainer(), output_var, **kwargs, shape=read_tensor_adaptor.Shape())

            # modify the write tensor data
            write_tensor_adaptor.data = read_tensor_adaptor.data * 2

            write_tensor_adaptor.StoreData()
            for entity in container:
                self.__CheckValues(getter_method(entity, input_var) * 2, getter_method(entity, output_var))

            # modify the read tensor data
            read_tensor_adaptor.data *= 6

            read_tensor_adaptor.StoreData()
            for entity in self.model_part.Nodes:
                self.__CheckValues(getter_method(entity, input_var), getter_method(entity, output_var) * 3)

    # def test_Dummy(self):
    #     import numpy

    #     a = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE, step_index=0, shape=[], is_historical=)
    # #     a.CollectData()

    #     b = VariableTensorAdaptor(a.GetContainer(), Kratos.DENSITY,)

    #     b.data = a.data

    #     # a_d = a.data
    #     # aa_d = a.ViewData()

    #     b.StoreData()
    #     for node in self.model_part.Nodes:
    #         print(node.GetSolutionStepValue(Kratos.DENSITY))

        # print(type(a_d))
        # print(type(aa_d))


    # def test_VariableTensorAdaptor(self):
        # nodal_tensor_adaptor = Kratos.TensorAdaptors.VariableTensorAdaptor()
        # nodal_tensor_adaptor.CollectData()
        # print(dir(Kratos.TensorAdaptors))
        pass

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
