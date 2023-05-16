import math

import KratosMultiphysics as Kratos
from KratosMultiphysics import IsDistributedRun
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.testing.utilities import ReadModelPart

class TestContainerExpressionUtils(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(KratosOA.HELMHOLTZ_VAR_DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        with kratos_unittest.WorkFolderScope(".", __file__, True):
            ReadModelPart("model_part_utils_test/quads", cls.model_part)

        for node in cls.model_part.Nodes:
            id = node.Id
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetSolutionStepValue(Kratos.PRESSURE, id+3)
            node.SetSolutionStepValue(Kratos.DENSITY, id+4)

    def test_ContainerVariableDataNormInf(self):
        a = Kratos.ContainerExpression.HistoricalExpression(self.model_part)

        a.Read(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerExpressionUtils.NormInf(a), 28)

        a.Read(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerExpressionUtils.NormInf(a), 30)

    def test_ContainerVariableDataNormL2(self):
        a = Kratos.ContainerExpression.HistoricalExpression(self.model_part)

        a.Read(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerExpressionUtils.NormL2(a), 87.74964387392122)

        a.Read(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerExpressionUtils.NormL2(a), 160.0781059358212)

    def test_ContainerVariableDataEntityMaxNormL2(self):
        a = Kratos.ContainerExpression.HistoricalExpression(self.model_part)

        a.Read(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerExpressionUtils.EntityMaxNormL2(a), 28)

        a.Read(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerExpressionUtils.EntityMaxNormL2(a), math.sqrt(28**2 + 29**2 + 30**2))

    def test_ContainerVariableDataInnerProduct(self):
        a = Kratos.ContainerExpression.HistoricalExpression(self.model_part)
        b = Kratos.ContainerExpression.HistoricalExpression(self.model_part)

        a.Read(Kratos.PRESSURE)
        b.Read(Kratos.DENSITY)

        self.assertEqual(KratosOA.ContainerExpressionUtils.InnerProduct(a, b), 8100.0)

    def test_ComputeNumberOfNeighbourConditions(self):
        neighbour_conditions = Kratos.ContainerExpression.NodalNonHistoricalExpression(self.model_part)
        KratosOA.ContainerExpressionUtils.ComputeNumberOfNeighbourConditions(neighbour_conditions)
        neighbour_conditions.Evaluate(Kratos.DENSITY)

        neighbour_map = {
            2.0: [1, 2, 4, 5, 7, 9, 10, 12, 15, 17, 20, 21, 22, 23, 24, 25],
            0.0: [3, 6, 8, 11, 13, 14, 16, 18, 19]
        }

        for node in self.model_part.Nodes:
            self.assertTrue(node.Id in neighbour_map[int(node.GetValue(Kratos.DENSITY))])

    def test_ComputeNumberOfNeighbourElements(self):
        neighbour_elements = Kratos.ContainerExpression.NodalNonHistoricalExpression(self.model_part)
        KratosOA.ContainerExpressionUtils.ComputeNumberOfNeighbourElements(neighbour_elements)
        neighbour_elements.Evaluate(Kratos.DENSITY)

        neighbour_map = {
            1.0: [1, 9, 22, 25],
            2.0: [2, 4, 5, 7, 10, 12, 15, 17, 20, 21, 23, 24],
            4.0: [3, 6, 8, 11, 13, 14, 16, 18, 19]
        }

        for node in self.model_part.Nodes:
            self.assertTrue(node.Id in neighbour_map[int(node.GetValue(Kratos.DENSITY))])

    def test_MapConditionVariableToNodalVariable(self):
        communicator: Kratos.Communicator = self.model_part.GetCommunicator()

        for condition in self.model_part.Conditions:
            condition.SetValue(Kratos.VELOCITY, Kratos.Array3([condition.Id, condition.Id + 1, condition.Id + 3]))
            condition.SetValue(Kratos.PRESSURE, condition.Id + 4)

        condition_container = Kratos.ContainerExpression.ConditionNonHistoricalExpression(self.model_part)
        neighbour_conditions = Kratos.ContainerExpression.NodalNonHistoricalExpression(self.model_part)
        mapped_values = Kratos.ContainerExpression.NodalNonHistoricalExpression(self.model_part)

        KratosOA.ContainerExpressionUtils.ComputeNumberOfNeighbourConditions(neighbour_conditions)
        neighbour_conditions.Evaluate(Kratos.YOUNG_MODULUS)

        condition_container.Read(Kratos.VELOCITY)
        KratosOA.ContainerExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, condition_container, neighbour_conditions)
        mapped_values.Evaluate(Kratos.VELOCITY)

        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.ACCELERATION, self.model_part.Nodes)
        for condition in self.model_part.Conditions:
            for node in condition.GetGeometry():
                node[Kratos.ACCELERATION] += condition.GetValue(Kratos.VELOCITY) / node.GetValue(Kratos.YOUNG_MODULUS)

        communicator.AssembleNonHistoricalData(Kratos.ACCELERATION)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node[Kratos.ACCELERATION], node[Kratos.VELOCITY])

        condition_container.Read(Kratos.PRESSURE)
        KratosOA.ContainerExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, condition_container, neighbour_conditions)
        mapped_values.Evaluate(Kratos.PRESSURE)

        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.DENSITY, self.model_part.Nodes)
        for condition in self.model_part.Conditions:
            for node in condition.GetGeometry():
                node[Kratos.DENSITY] += condition.GetValue(Kratos.PRESSURE) / node.GetValue(Kratos.YOUNG_MODULUS)

        communicator.AssembleNonHistoricalData(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertAlmostEqual(node[Kratos.DENSITY], node[Kratos.PRESSURE], 9)

    def test_MapElementVariableToNodalVariable(self):
        communicator: Kratos.Communicator = self.model_part.GetCommunicator()

        for element in self.model_part.Elements:
            element.SetValue(Kratos.VELOCITY, Kratos.Array3([element.Id, element.Id + 1, element.Id + 3]))
            element.SetValue(Kratos.PRESSURE, element.Id + 4)

        element_container = Kratos.ContainerExpression.ElementNonHistoricalExpression(self.model_part)
        neighbour_conditions = Kratos.ContainerExpression.NodalNonHistoricalExpression(self.model_part)
        mapped_values = Kratos.ContainerExpression.NodalNonHistoricalExpression(self.model_part)

        KratosOA.ContainerExpressionUtils.ComputeNumberOfNeighbourElements(neighbour_conditions)
        neighbour_conditions.Evaluate(Kratos.YOUNG_MODULUS)

        element_container.Read(Kratos.VELOCITY)
        KratosOA.ContainerExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, element_container, neighbour_conditions)
        mapped_values.Evaluate(Kratos.VELOCITY)

        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.ACCELERATION, self.model_part.Nodes)
        for element in self.model_part.Elements:
            for node in element.GetGeometry():
                node[Kratos.ACCELERATION] += element.GetValue(Kratos.VELOCITY) / node.GetValue(Kratos.YOUNG_MODULUS)

        communicator.AssembleNonHistoricalData(Kratos.ACCELERATION)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node[Kratos.ACCELERATION], node[Kratos.VELOCITY])

        element_container.Read(Kratos.PRESSURE)
        KratosOA.ContainerExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, element_container, neighbour_conditions)
        mapped_values.Evaluate(Kratos.PRESSURE)

        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.DENSITY, self.model_part.Nodes)
        for element in self.model_part.Elements:
            for node in element.GetGeometry():
                node[Kratos.DENSITY] += element.GetValue(Kratos.PRESSURE) / node.GetValue(Kratos.YOUNG_MODULUS)

        communicator.AssembleNonHistoricalData(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertEqual(node[Kratos.DENSITY], node[Kratos.PRESSURE])

    def test_MapNodalVariableToConditionVariable(self):
        for node in self.model_part.Nodes:
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([node.Id, node.Id + 1, node.Id + 3]))
            node.SetValue(Kratos.PRESSURE, node.Id + 4)

        nodal_container = Kratos.ContainerExpression.NodalNonHistoricalExpression(self.model_part)
        mapped_value = Kratos.ContainerExpression.ConditionNonHistoricalExpression(self.model_part)

        nodal_container.Read(Kratos.VELOCITY)
        KratosOA.ContainerExpressionUtils.MapNodalVariableToContainerVariable(mapped_value, nodal_container)
        mapped_value.Evaluate(Kratos.ACCELERATION)

        for condition in self.model_part.Conditions:
            v = Kratos.Array3([0, 0, 0])
            for node in condition.GetGeometry():
                v += node.GetValue(Kratos.VELOCITY)
            self.assertVectorAlmostEqual(v / 2.0, condition.GetValue(Kratos.ACCELERATION))

        nodal_container.Read(Kratos.PRESSURE)
        KratosOA.ContainerExpressionUtils.MapNodalVariableToContainerVariable(mapped_value, nodal_container)
        mapped_value.Evaluate(Kratos.DENSITY)

        for condition in self.model_part.Conditions:
            v = 0.0
            for node in condition.GetGeometry():
                v += node.GetValue(Kratos.PRESSURE)
            self.assertEqual(v / 2.0, condition.GetValue(Kratos.DENSITY))

    def test_ProductWithEntityMatrix(self):
        if (IsDistributedRun()):
            self.skipTest("Skipping since ProductWithEntityMatrix does not support MPI yet.")

        number_of_nodes = self.model_part.NumberOfNodes()

        a = Kratos.ContainerExpression.HistoricalExpression(self.model_part)
        a.Read(Kratos.PRESSURE)

        m = Kratos.Matrix(number_of_nodes, number_of_nodes)
        for i in range(number_of_nodes):
            for j in range(number_of_nodes):
                m[i, j] = (i + 1) * (j + 1)

        b = Kratos.ContainerExpression.HistoricalExpression(self.model_part)
        KratosOA.ContainerExpressionUtils.ProductWithEntityMatrix(b, m, a)
        b.Evaluate(Kratos.DENSITY)

        for i, node_b in enumerate(b.GetContainer()):
            v = 0
            for j, node_a in enumerate(a.GetContainer()):
                v += m[i, j] * node_a.GetSolutionStepValue(Kratos.PRESSURE)
            self.assertEqual(v, node_b.GetSolutionStepValue(Kratos.DENSITY))

    def test_ProductWithEntityMatrixSparse(self):
        if (IsDistributedRun()):
            self.skipTest("Skipping since ProductWithEntityMatrix does not support MPI yet.")

        number_of_nodes = self.model_part.NumberOfNodes()

        a = Kratos.ContainerExpression.HistoricalExpression(self.model_part)
        a.Read(Kratos.PRESSURE)

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

        dense_b = Kratos.ContainerExpression.HistoricalExpression(self.model_part)
        KratosOA.ContainerExpressionUtils.ProductWithEntityMatrix(dense_b, dense_m, a)

        sparse_b = Kratos.ContainerExpression.HistoricalExpression(self.model_part)
        KratosOA.ContainerExpressionUtils.ProductWithEntityMatrix(sparse_b, sparse_m, a)

        self.assertEqual(KratosOA.ContainerExpressionUtils.InnerProduct(dense_b - sparse_b, dense_b - sparse_b), 0)

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
        KratosOA.ContainerExpressionUtils.Transpose(transpose_dense_m, dense_m)
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
        KratosOA.ContainerExpressionUtils.Transpose(transpose_sparse_m, sparse_m)
        for i in range(dense_m.Size1()):
            for j in range(dense_m.Size2()):
                self.assertEqual(transpose_sparse_m[j, i], dense_m[i, j])

    def test_ComputeVariableDataHolderProductWithEntityMatrix(self):
        for node in self.model_part.Nodes:
            node.SetValue(Kratos.PRESSURE, node.Id + 1)
            node.SetSolutionStepValue(KratosOA.HELMHOLTZ_VAR_DENSITY, node.Id + 1)

        nodal_values = Kratos.ContainerExpression.NodalNonHistoricalExpression(self.model_part)
        nodal_values.Read(Kratos.PRESSURE)

        output_values = Kratos.ContainerExpression.NodalNonHistoricalExpression(self.model_part)
        KratosOA.ContainerExpressionUtils.ComputeNodalVariableProductWithEntityMatrix(output_values, nodal_values, KratosOA.HELMHOLTZ_MASS_MATRIX, self.model_part.Elements)

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

        analytical_values = Kratos.ContainerExpression.NodalNonHistoricalExpression(self.model_part)
        analytical_values.Read(Kratos.DENSITY)

        self.assertEqual(KratosOA.ContainerExpressionUtils.NormL2(analytical_values - output_values), 0.0)

    def test_Scopes(self):
        a = Kratos.ContainerExpression.HistoricalExpression(self.model_part)
        b = Kratos.ContainerExpression.HistoricalExpression(self.model_part)

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.PRESSURE)

        def func(a):
            c = a * 2
            d = c * 3
            return d

        c = func(a)
        d = func(b)

        e = Kratos.ContainerExpression.NodalNonHistoricalExpression(c)
        e.Evaluate(Kratos.ACCELERATION)
        f = Kratos.ContainerExpression.NodalNonHistoricalExpression(d)
        f.Evaluate(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.VELOCITY) * 6, node.GetValue(Kratos.ACCELERATION))
            self.assertEqual(node.GetSolutionStepValue(Kratos.PRESSURE) * 6, node.GetValue(Kratos.DENSITY))

    def test_CollectiveExpressionsNormInf(self):
        a = Kratos.ContainerExpression.HistoricalExpression(self.model_part)
        b = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.PRESSURE)

        collective_1 = KratosOA.ContainerExpression.CollectiveExpressions([a, b])
        self.assertEqual(KratosOA.ContainerExpressionUtils.NormInf(collective_1), max(KratosOA.ContainerExpressionUtils.NormInf(a), KratosOA.ContainerExpressionUtils.NormInf(b)))

    def test_CollectiveExpressionsNormL2(self):
        a = Kratos.ContainerExpression.HistoricalExpression(self.model_part)
        b = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.PRESSURE)

        collective_1 = KratosOA.ContainerExpression.CollectiveExpressions([a, b])
        self.assertEqual(KratosOA.ContainerExpressionUtils.NormL2(collective_1), math.sqrt(KratosOA.ContainerExpressionUtils.NormL2(a)**2 + KratosOA.ContainerExpressionUtils.NormL2(b)**2))


    def test_CollectiveExpressionsInnerProduct(self):
        a = Kratos.ContainerExpression.HistoricalExpression(self.model_part)
        b = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)

        collective_1 = KratosOA.ContainerExpression.CollectiveExpressions([a, b])

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.PRESSURE)

        self.assertEqual(KratosOA.ContainerExpressionUtils.InnerProduct(collective_1, collective_1), KratosOA.ContainerExpressionUtils.InnerProduct(a, a) + KratosOA.ContainerExpressionUtils.InnerProduct(b, b))

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()