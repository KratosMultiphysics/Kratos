
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

    def setUp(self):
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.Id + 1)
            node.SetSolutionStepValue(Kratos.VELOCITY, [node.Id * 2 + 1, node.Id * 3 + 1, node.Id * 4 + 1])
            node.SetSolutionStepValue(Kratos.EXTERNAL_FORCES_VECTOR, Kratos.Vector([node.Id * 2 + 1, node.Id * 3 + 1, node.Id * 4 + 1, node.Id * 2 + 1, node.Id * 3 + 1]))
            node.SetSolutionStepValue(Kratos.NORMAL_SHAPE_DERIVATIVE, Kratos.Matrix([[node.Id * 2 + 1, node.Id * 3 + 1], [node.Id * 4 + 1, node.Id * 2 + 1], [node.Id * 3 + 1, node.Id * 4 + 1]]))

    def test_NodeHistoricalVariableTensorAdaptor(self):
        for input_var, output_var in zip(self.input_list_of_variables, self.output_list_of_variables):
            read_tensor_adaptor = Kratos.TensorAdaptors.NodeTensorAdaptors.NodeHistoricalVariableTensorAdaptor(self.model_part.Nodes, input_var, 0)
            write_tensor_adaptor = Kratos.TensorAdaptors.NodeTensorAdaptors.NodeHistoricalVariableTensorAdaptor(read_tensor_adaptor.GetContainer(), output_var, read_tensor_adaptor.Shape(), 1)

            read_tensor_adaptor.CollectData()

            # modify the write tensor data
            write_tensor_adaptor.data = read_tensor_adaptor.data * 2

            write_tensor_adaptor.StoreData()
            for node in self.model_part.Nodes:
                self.__CheckValues(node.GetSolutionStepValue(input_var) * 2, node.GetSolutionStepValue(output_var, 1))

            # modify the read tensor data
            read_tensor_adaptor.data *= 6

            read_tensor_adaptor.StoreData()
            for node in self.model_part.Nodes:
                self.__CheckValues(node.GetSolutionStepValue(input_var), node.GetSolutionStepValue(output_var, 1) * 3)

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
