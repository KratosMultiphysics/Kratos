import math

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest


class TestContainerVariableDataHolderUtils(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        cls.number_of_nodes = 10
        for id in range(1, cls.number_of_nodes + 1):
            node = cls.model_part.CreateNewNode(id, id, id+1, id+2)
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetSolutionStepValue(Kratos.PRESSURE, id+3)
            node.SetSolutionStepValue(Kratos.DENSITY, id+4)

    def test_NormInf(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormInf(a), 13)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormInf(a), 15)

    def test_EntityMaxNormL2(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.EntityMaxNormL2(a), 13)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.EntityMaxNormL2(a), math.sqrt(15**2 + 14**2 + 13**2))

    def test_InnerProduct(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        b.ReadDataFromContainerVariable(Kratos.DENSITY)

        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.InnerProduct(a, b), 890)

    def test_ProductWithEntityMatrix(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        a.ReadDataFromContainerVariable(Kratos.PRESSURE)

        m = Kratos.Matrix(self.number_of_nodes, self.number_of_nodes)
        for i in range(self.number_of_nodes):
            for j in range(self.number_of_nodes):
                m[i, j] = (i + 1) * (j + 1)

        b = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        KratosOA.ContainerVariableDataHolderUtils.ProductWithEntityMatrix(b, m, a)
        b.AssignDataToContainerVariable(Kratos.DENSITY)

        for i, node_b in enumerate(b.GetContainer()):
            v = 0
            for j, node_a in enumerate(a.GetContainer()):
                v += m[i, j] * node_a.GetSolutionStepValue(Kratos.PRESSURE)
            self.assertEqual(v, node_b.GetSolutionStepValue(Kratos.DENSITY))

    def test_ProductWithEntityMatrixSparse(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        a.ReadDataFromContainerVariable(Kratos.PRESSURE)

        # first build the normal matrix
        dense_m = Kratos.Matrix(
            [
                [5, 0, 0, 2, 0, 0, 0, 0, 0, 2],
                [0, 0, 3, 2, 3, 0, 0, 0, 0, 0],
                [0, 2, 0, 2, 0, 2, 0, 6, 0, 0],
                [0, 0, 0, 2, 0, 4, 0, 5, 0, 0],
                [0, 0, 2, 0, 0, 2, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 3, 3, 0],
                [0, 2, 0, 0, 3, 0, 4, 0, 0, 0],
                [3, 0, 0, 2, 0, 0, 0, 2, 0, 2],
                [9, 0, 0, 0, 0, 7, 0, 7, 0, 5],
                [0, 0, 7, 2, 0, 0, 0, 0, 0, 2]
            ]
        )

        # now add values to the sparse matrix
        sparse_m = Kratos.CompressedMatrix(self.number_of_nodes, self.number_of_nodes)
        for i in range(self.number_of_nodes):
            for j in range(self.number_of_nodes):
                if dense_m[i, j] != 0.0:
                    sparse_m[i, j] = dense_m[i, j]

        dense_b = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        KratosOA.ContainerVariableDataHolderUtils.ProductWithEntityMatrix(dense_b, dense_m, a)

        sparse_b = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        KratosOA.ContainerVariableDataHolderUtils.ProductWithEntityMatrix(sparse_b, sparse_m, a)

        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.InnerProduct(dense_b - sparse_b, dense_b - sparse_b), 0)

    def test_Transpose(self):
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
        KratosOA.ContainerVariableDataHolderUtils.Transpose(transpose_dense_m, dense_m)
        for i in range(dense_m.Size1()):
            for j in range(dense_m.Size2()):
                self.assertEqual(transpose_dense_m[j, i], dense_m[i, j])

        # now add values to the sparse matrix
        sparse_m = Kratos.CompressedMatrix(self.number_of_nodes, self.number_of_nodes)
        for i in range(dense_m.Size1()):
            for j in range(dense_m.Size2()):
                if dense_m[i, j] != 0.0:
                    sparse_m[i, j] = dense_m[i, j]

        transpose_sparse_m = Kratos.CompressedMatrix()
        KratosOA.ContainerVariableDataHolderUtils.Transpose(transpose_sparse_m, sparse_m)
        for i in range(dense_m.Size1()):
            for j in range(dense_m.Size2()):
                self.assertEqual(transpose_sparse_m[j, i], dense_m[i, j])

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()