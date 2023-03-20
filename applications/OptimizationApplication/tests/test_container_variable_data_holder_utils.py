import math

import KratosMultiphysics as Kratos
from KratosMultiphysics import IsDistributedRun
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.testing.utilities import ReadModelPart

class TestContainerVariableDataUtils(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(KratosOA.HELMHOLTZ_VAR_DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        ReadModelPart("model_part_utils_test/quads", cls.model_part)

        for node in cls.model_part.Nodes:
            id = node.Id
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetSolutionStepValue(Kratos.PRESSURE, id+3)
            node.SetSolutionStepValue(Kratos.DENSITY, id+4)

    def test_ContainerVariableDataNormInf(self):
        a = KratosOA.HistoricalVariableData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerVariableDataUtils.NormInf(a), 28)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerVariableDataUtils.NormInf(a), 30)

    def test_ContainerVariableDataNormL2(self):
        a = KratosOA.HistoricalVariableData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerVariableDataUtils.NormL2(a), 87.74964387392122)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerVariableDataUtils.NormL2(a), 160.0781059358212)

    def test_ContainerVariableDataEntityMaxNormL2(self):
        a = KratosOA.HistoricalVariableData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerVariableDataUtils.EntityMaxNormL2(a), 28)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerVariableDataUtils.EntityMaxNormL2(a), math.sqrt(28**2 + 29**2 + 30**2))

    def test_ContainerVariableDataInnerProduct(self):
        a = KratosOA.HistoricalVariableData(self.model_part)
        b = KratosOA.HistoricalVariableData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        b.ReadDataFromContainerVariable(Kratos.DENSITY)

        self.assertEqual(KratosOA.ContainerVariableDataUtils.InnerProduct(a, b), 8100.0)

    def test_ComputeNumberOfNeighbourConditions(self):
        neighbour_conditions = KratosOA.NodalNonHistoricalVariableData(self.model_part)
        KratosOA.ContainerVariableDataUtils.ComputeNumberOfNeighbourConditions(neighbour_conditions)
        neighbour_conditions.AssignDataToContainerVariable(Kratos.DENSITY)

        neighbour_map = {
            2.0: [1, 2, 4, 5, 7, 9, 10, 12, 15, 17, 20, 21, 22, 23, 24, 25],
            0.0: [3, 6, 8, 11, 13, 14, 16, 18, 19]
        }

        for node in self.model_part.Nodes:
            self.assertTrue(node.Id in neighbour_map[int(node.GetValue(Kratos.DENSITY))])

    def test_ComputeNumberOfNeighbourElements(self):
        neighbour_elements = KratosOA.NodalNonHistoricalVariableData(self.model_part)
        KratosOA.ContainerVariableDataUtils.ComputeNumberOfNeighbourElements(neighbour_elements)
        neighbour_elements.AssignDataToContainerVariable(Kratos.DENSITY)

        neighbour_map = {
            1.0: [1, 9, 22, 25],
            2.0: [2, 4, 5, 7, 10, 12, 15, 17, 20, 21, 23, 24],
            4.0: [3, 6, 8, 11, 13, 14, 16, 18, 19]
        }

        for node in self.model_part.Nodes:
            self.assertTrue(node.Id in neighbour_map[int(node.GetValue(Kratos.DENSITY))])

    def test_MapConditionContainerVariableDataToNodalVariableDataHolder(self):
        communicator: Kratos.Communicator = self.model_part.GetCommunicator()

        for condition in self.model_part.Conditions:
            condition.SetValue(Kratos.VELOCITY, Kratos.Array3([condition.Id, condition.Id + 1, condition.Id + 3]))
            condition.SetValue(Kratos.PRESSURE, condition.Id + 4)

        condition_container = KratosOA.ConditionNonHistoricalVariableData(self.model_part)
        neighbour_conditions = KratosOA.NodalNonHistoricalVariableData(self.model_part)
        mapped_values = KratosOA.NodalNonHistoricalVariableData(self.model_part)

        KratosOA.ContainerVariableDataUtils.ComputeNumberOfNeighbourConditions(neighbour_conditions)
        neighbour_conditions.AssignDataToContainerVariable(Kratos.YOUNG_MODULUS)

        condition_container.ReadDataFromContainerVariable(Kratos.VELOCITY)
        KratosOA.ContainerVariableDataUtils.MapContainerVariableDataToNodalVariableData(mapped_values, condition_container, neighbour_conditions)
        mapped_values.AssignDataToContainerVariable(Kratos.VELOCITY)

        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.ACCELERATION, self.model_part.Nodes)
        for condition in self.model_part.Conditions:
            for node in condition.GetGeometry():
                node[Kratos.ACCELERATION] += condition.GetValue(Kratos.VELOCITY) / node.GetValue(Kratos.YOUNG_MODULUS)

        communicator.AssembleNonHistoricalData(Kratos.ACCELERATION)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node[Kratos.ACCELERATION], node[Kratos.VELOCITY])

        condition_container.ReadDataFromContainerVariable(Kratos.PRESSURE)
        KratosOA.ContainerVariableDataUtils.MapContainerVariableDataToNodalVariableData(mapped_values, condition_container, neighbour_conditions)
        mapped_values.AssignDataToContainerVariable(Kratos.PRESSURE)

        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.DENSITY, self.model_part.Nodes)
        for condition in self.model_part.Conditions:
            for node in condition.GetGeometry():
                node[Kratos.DENSITY] += condition.GetValue(Kratos.PRESSURE) / node.GetValue(Kratos.YOUNG_MODULUS)

        communicator.AssembleNonHistoricalData(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertAlmostEqual(node[Kratos.DENSITY], node[Kratos.PRESSURE], 9)

    def test_MapElementContainerVariableDataToNodalVariableDataHolder(self):
        communicator: Kratos.Communicator = self.model_part.GetCommunicator()

        for element in self.model_part.Elements:
            element.SetValue(Kratos.VELOCITY, Kratos.Array3([element.Id, element.Id + 1, element.Id + 3]))
            element.SetValue(Kratos.PRESSURE, element.Id + 4)

        element_container = KratosOA.ElementNonHistoricalVariableData(self.model_part)
        neighbour_conditions = KratosOA.NodalNonHistoricalVariableData(self.model_part)
        mapped_values = KratosOA.NodalNonHistoricalVariableData(self.model_part)

        KratosOA.ContainerVariableDataUtils.ComputeNumberOfNeighbourElements(neighbour_conditions)
        neighbour_conditions.AssignDataToContainerVariable(Kratos.YOUNG_MODULUS)

        element_container.ReadDataFromContainerVariable(Kratos.VELOCITY)
        KratosOA.ContainerVariableDataUtils.MapContainerVariableDataToNodalVariableData(mapped_values, element_container, neighbour_conditions)
        mapped_values.AssignDataToContainerVariable(Kratos.VELOCITY)

        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.ACCELERATION, self.model_part.Nodes)
        for element in self.model_part.Elements:
            for node in element.GetGeometry():
                node[Kratos.ACCELERATION] += element.GetValue(Kratos.VELOCITY) / node.GetValue(Kratos.YOUNG_MODULUS)

        communicator.AssembleNonHistoricalData(Kratos.ACCELERATION)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node[Kratos.ACCELERATION], node[Kratos.VELOCITY])

        element_container.ReadDataFromContainerVariable(Kratos.PRESSURE)
        KratosOA.ContainerVariableDataUtils.MapContainerVariableDataToNodalVariableData(mapped_values, element_container, neighbour_conditions)
        mapped_values.AssignDataToContainerVariable(Kratos.PRESSURE)

        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.DENSITY, self.model_part.Nodes)
        for element in self.model_part.Elements:
            for node in element.GetGeometry():
                node[Kratos.DENSITY] += element.GetValue(Kratos.PRESSURE) / node.GetValue(Kratos.YOUNG_MODULUS)

        communicator.AssembleNonHistoricalData(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertEqual(node[Kratos.DENSITY], node[Kratos.PRESSURE])

    def test_MapNodalVariableDataHolderToConditionContainerVariableDataHolder(self):
        for node in self.model_part.Nodes:
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([node.Id, node.Id + 1, node.Id + 3]))
            node.SetValue(Kratos.PRESSURE, node.Id + 4)

        nodal_container = KratosOA.NodalNonHistoricalVariableData(self.model_part)
        mapped_value = KratosOA.ConditionNonHistoricalVariableData(self.model_part)

        nodal_container.ReadDataFromContainerVariable(Kratos.VELOCITY)
        KratosOA.ContainerVariableDataUtils.MapNodalVariableDataToContainerVariableData(mapped_value, nodal_container)
        mapped_value.AssignDataToContainerVariable(Kratos.ACCELERATION)

        for condition in self.model_part.Conditions:
            v = Kratos.Array3([0, 0, 0])
            for node in condition.GetGeometry():
                v += node.GetValue(Kratos.VELOCITY)
            self.assertVectorAlmostEqual(v / 2.0, condition.GetValue(Kratos.ACCELERATION))

        nodal_container.ReadDataFromContainerVariable(Kratos.PRESSURE)
        KratosOA.ContainerVariableDataUtils.MapNodalVariableDataToContainerVariableData(mapped_value, nodal_container)
        mapped_value.AssignDataToContainerVariable(Kratos.DENSITY)

        for condition in self.model_part.Conditions:
            v = 0.0
            for node in condition.GetGeometry():
                v += node.GetValue(Kratos.PRESSURE)
            self.assertEqual(v / 2.0, condition.GetValue(Kratos.DENSITY))

    def test_ProductWithEntityMatrix(self):
        if (IsDistributedRun()):
            self.skipTest("Skipping since ProductWithEntityMatrix does not support MPI yet.")

        number_of_nodes = self.model_part.NumberOfNodes()

        a = KratosOA.HistoricalVariableData(self.model_part)
        a.ReadDataFromContainerVariable(Kratos.PRESSURE)

        m = Kratos.Matrix(number_of_nodes, number_of_nodes)
        for i in range(number_of_nodes):
            for j in range(number_of_nodes):
                m[i, j] = (i + 1) * (j + 1)

        b = KratosOA.HistoricalVariableData(self.model_part)
        KratosOA.ContainerVariableDataUtils.ProductWithEntityMatrix(b, m, a)
        b.AssignDataToContainerVariable(Kratos.DENSITY)

        for i, node_b in enumerate(b.GetContainer()):
            v = 0
            for j, node_a in enumerate(a.GetContainer()):
                v += m[i, j] * node_a.GetSolutionStepValue(Kratos.PRESSURE)
            self.assertEqual(v, node_b.GetSolutionStepValue(Kratos.DENSITY))

    def test_ProductWithEntityMatrixSparse(self):
        if (IsDistributedRun()):
            self.skipTest("Skipping since ProductWithEntityMatrix does not support MPI yet.")

        number_of_nodes = self.model_part.NumberOfNodes()

        a = KratosOA.HistoricalVariableData(self.model_part)
        a.ReadDataFromContainerVariable(Kratos.PRESSURE)

        dense_m = Kratos.Matrix(number_of_nodes, number_of_nodes)
        for i in range(number_of_nodes):
            for j in range(number_of_nodes):
                dense_m[i, j] = int((i + 1) * (j + 1) % 10)

        # now add values to the sparse matrix
        sparse_m = Kratos.CompressedMatrix(number_of_nodes, number_of_nodes)
        for i in range(number_of_nodes):
            for j in range(number_of_nodes):
                if dense_m[i, j] != 0.0:
                    sparse_m[i, j] = dense_m[i, j]

        dense_b = KratosOA.HistoricalVariableData(self.model_part)
        KratosOA.ContainerVariableDataUtils.ProductWithEntityMatrix(dense_b, dense_m, a)

        sparse_b = KratosOA.HistoricalVariableData(self.model_part)
        KratosOA.ContainerVariableDataUtils.ProductWithEntityMatrix(sparse_b, sparse_m, a)

        self.assertEqual(KratosOA.ContainerVariableDataUtils.InnerProduct(dense_b - sparse_b, dense_b - sparse_b), 0)

    def test_Transpose(self):
        matrix_size = 10

        # first build the normal matrix
        dense_m = Kratos.Matrix(
            [
                [5, 0, 0, 2, 0, 0, 0, 0, 0],
                [0, 0, 3, 2, 3, 0, 0, 0, 0],
                [0, 2, 0, 2, 0, 2, 0, 6, 0],
                [0, 0, 0, 2, 0, 4, 0, 5, 0],
                [0, 0, 2, 0, 0, 2, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 3, 3],
                [0, 2, 0, 0, 3, 0, 4, 0, 0],
                [3, 0, 0, 2, 0, 0, 0, 2, 0],
                [9, 0, 0, 0, 0, 7, 0, 7, 0],
                [0, 0, 7, 2, 0, 0, 0, 0, 0]
            ]
        )
        transpose_dense_m = Kratos.Matrix()
        KratosOA.ContainerVariableDataUtils.Transpose(transpose_dense_m, dense_m)
        for i in range(dense_m.Size1()):
            for j in range(dense_m.Size2()):
                self.assertEqual(transpose_dense_m[j, i], dense_m[i, j])

        # now add values to the sparse matrix
        sparse_m = Kratos.CompressedMatrix(matrix_size, matrix_size)
        for i in range(dense_m.Size1()):
            for j in range(dense_m.Size2()):
                if dense_m[i, j] != 0.0:
                    sparse_m[i, j] = dense_m[i, j]

        transpose_sparse_m = Kratos.CompressedMatrix()
        KratosOA.ContainerVariableDataUtils.Transpose(transpose_sparse_m, sparse_m)
        for i in range(dense_m.Size1()):
            for j in range(dense_m.Size2()):
                self.assertEqual(transpose_sparse_m[j, i], dense_m[i, j])

    def test_ComputeVariableDataHolderProductWithEntityMatrix(self):
        for node in self.model_part.Nodes:
            node.SetValue(Kratos.PRESSURE, node.Id + 1)
            node.SetSolutionStepValue(KratosOA.HELMHOLTZ_VAR_DENSITY, node.Id + 1)

        nodal_values = KratosOA.NodalNonHistoricalVariableData(self.model_part)
        nodal_values.ReadDataFromContainerVariable(Kratos.PRESSURE)

        output_values = KratosOA.NodalNonHistoricalVariableData(self.model_part)
        KratosOA.ContainerVariableDataUtils.ComputeVariableDataProductWithEntityMatrix(output_values, nodal_values, KratosOA.HELMHOLTZ_MASS_MATRIX, self.model_part.Elements)

        # analytical calculation
        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.DENSITY, self.model_part.Nodes)
        element: Kratos.Element
        for element in self.model_part.Elements:
            v_in = Kratos.Vector(4)
            for i, node in enumerate(element.GetGeometry()):
                v_in[i] = node.GetSolutionStepValue(KratosOA.HELMHOLTZ_VAR_DENSITY)

            m = element.Calculate(KratosOA.HELMHOLTZ_MASS_MATRIX, self.model_part.ProcessInfo)
            v_out = m * v_in

            for i, node in enumerate(element.GetGeometry()):
                node.SetValue(Kratos.DENSITY, node.GetValue(Kratos.DENSITY) + v_out[i])

        self.model_part.GetCommunicator().AssembleNonHistoricalData(Kratos.DENSITY)

        analytical_values = KratosOA.NodalNonHistoricalVariableData(self.model_part)
        analytical_values.ReadDataFromContainerVariable(Kratos.DENSITY)

        self.assertEqual(KratosOA.ContainerVariableDataUtils.NormL2(analytical_values - output_values), 0.0)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()