import math

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
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

        number_of_nodes = 10
        for id in range(1, number_of_nodes + 1):
            node = cls.model_part.CreateNewNode(id, id, id+1, id+2)
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetSolutionStepValue(Kratos.PRESSURE, id+3)
            node.SetSolutionStepValue(Kratos.DENSITY, id+4)

    @staticmethod
    def __GetNodeIndex(position, model_part):
        for node in model_part.Nodes:
            distance = (position[0] - node.X)**2 + (position[1] - node.Y)**2 + (position[2] - node.Z)**2
            if distance < 1e-8:
                return node.Id

        node = model_part.CreateNewNode(model_part.NumberOfNodes(), position[0], position[1], position[2])
        return node.Id

    def test_ContainerVariableDataHolderNormInf(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormInf(a), 13)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormInf(a), 15)

    def test_ContainerVariableDataHolderNormL2(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)

        l2_norm = 0.0
        for node in self.model_part.Nodes:
            l2_norm += node.GetSolutionStepValue(Kratos.PRESSURE) ** 2
        l2_norm = math.sqrt(l2_norm)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormL2(a), l2_norm)

        l2_norm = 0.0
        for node in self.model_part.Nodes:
            v = node.GetSolutionStepValue(Kratos.VELOCITY)
            l2_norm += v[0] ** 2 + v[1] ** 2 + v[2] ** 2
        l2_norm = math.sqrt(l2_norm)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormL2(a), l2_norm)

    def test_ContainerVariableDataHolderEntityMaxNormL2(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.EntityMaxNormL2(a), 13)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.EntityMaxNormL2(a), math.sqrt(15**2 + 14**2 + 13**2))

    def test_ContainerVariableDataHolderInnerProduct(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        b.ReadDataFromContainerVariable(Kratos.DENSITY)

        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.InnerProduct(a, b), 890)

    def test_CollectiveVariableDataHolderNormInf(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveVariableDataHolder([a, b])
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormInf(collective_1), max(KratosOA.ContainerVariableDataHolderUtils.NormInf(a), KratosOA.ContainerVariableDataHolderUtils.NormInf(b)))

    def test_CollectiveVariableDataHolderNormL2(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveVariableDataHolder([a, b])
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormL2(collective_1), math.sqrt(KratosOA.ContainerVariableDataHolderUtils.NormL2(a)**2 + KratosOA.ContainerVariableDataHolderUtils.NormL2(b)**2))


    def test_CollectiveVariableDataHolderInnerProduct(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        collective_1 = KratosOA.CollectiveVariableDataHolder([a, b])

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.InnerProduct(collective_1, collective_1), KratosOA.ContainerVariableDataHolderUtils.InnerProduct(a, a) + KratosOA.ContainerVariableDataHolderUtils.InnerProduct(b, b))

    def test_ProductWithEntityMatrix(self):
        number_of_nodes = self.model_part.NumberOfNodes()

        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        a.ReadDataFromContainerVariable(Kratos.PRESSURE)

        m = Kratos.Matrix(number_of_nodes, number_of_nodes)
        for i in range(number_of_nodes):
            for j in range(number_of_nodes):
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
        number_of_nodes = self.model_part.NumberOfNodes()

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
        sparse_m = Kratos.CompressedMatrix(number_of_nodes, number_of_nodes)
        for i in range(number_of_nodes):
            for j in range(number_of_nodes):
                if dense_m[i, j] != 0.0:
                    sparse_m[i, j] = dense_m[i, j]

        dense_b = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        KratosOA.ContainerVariableDataHolderUtils.ProductWithEntityMatrix(dense_b, dense_m, a)

        sparse_b = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        KratosOA.ContainerVariableDataHolderUtils.ProductWithEntityMatrix(sparse_b, sparse_m, a)

        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.InnerProduct(dense_b - sparse_b, dense_b - sparse_b), 0)

    def test_Transpose(self):
        number_of_nodes = self.model_part.NumberOfNodes()

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
        sparse_m = Kratos.CompressedMatrix(number_of_nodes, number_of_nodes)
        for i in range(dense_m.Size1()):
            for j in range(dense_m.Size2()):
                if dense_m[i, j] != 0.0:
                    sparse_m[i, j] = dense_m[i, j]

        transpose_sparse_m = Kratos.CompressedMatrix()
        KratosOA.ContainerVariableDataHolderUtils.Transpose(transpose_sparse_m, sparse_m)
        for i in range(dense_m.Size1()):
            for j in range(dense_m.Size2()):
                self.assertEqual(transpose_sparse_m[j, i], dense_m[i, j])

    def test_ComputeNumberOfNeighbourConditions(self):
        def __GenerateEntity(center: Kratos.Array3, model_part: Kratos.ModelPart):
            model_part.CreateNewCondition("SurfaceCondition3D4N", model_part.NumberOfConditions() + 1, [
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([-0.5, -0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([0.5, -0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([0.5, 0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([-0.5, 0.5, 0.0]), model_part)
            ], model_part.GetProperties()[1])

        model = Kratos.Model()
        model_part = model.CreateModelPart("test1")
        model_part.CreateNewProperties(1)
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2

        for i in range(9):
            __GenerateEntity(Kratos.Array3([i // 3, i % 3, 0]), model_part)

        dummy_condition_container = KratosOA.ConditionContainerVariableDataHolder(model_part)
        neighbour_elements = KratosOA.NodalContainerVariableDataHolder(model_part)
        KratosOA.ContainerVariableDataHolderUtils.ComputeNumberOfNeighbourEntities(neighbour_elements, dummy_condition_container)
        neighbour_elements.AssignDataToContainerVariable(Kratos.DENSITY)

        neighbour_map = {
            1: [1, 8, 13, 16],
            2: [2, 9, 4, 6, 14, 15, 7, 12],
            4: [3, 10, 11, 5]
        }
        for node in model_part.Nodes:
            self.assertTrue(node.Id + 1 in neighbour_map[int(node.GetValue(Kratos.DENSITY))])

    def test_ComputeNumberOfNeighbourConditions(self):
        def __GenerateEntity(center: Kratos.Array3, model_part: Kratos.ModelPart):
            model_part.CreateNewElement("Element2D4N", model_part.NumberOfElements() + 1, [
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([-0.5, -0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([0.5, -0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([0.5, 0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([-0.5, 0.5, 0.0]), model_part)
            ], model_part.GetProperties()[1])

        model = Kratos.Model()
        model_part = model.CreateModelPart("test1")
        model_part.CreateNewProperties(1)
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2

        for i in range(9):
            __GenerateEntity(Kratos.Array3([i // 3, i % 3, 0]), model_part)

        dummy_condition_container = KratosOA.ElementContainerVariableDataHolder(model_part)
        neighbour_elements = KratosOA.NodalContainerVariableDataHolder(model_part)
        KratosOA.ContainerVariableDataHolderUtils.ComputeNumberOfNeighbourEntities(neighbour_elements, dummy_condition_container)
        neighbour_elements.AssignDataToContainerVariable(Kratos.DENSITY)

        neighbour_map = {
            1: [1, 8, 13, 16],
            2: [2, 9, 4, 6, 14, 15, 7, 12],
            4: [3, 10, 11, 5]
        }
        for node in model_part.Nodes:
            self.assertTrue(node.Id + 1 in neighbour_map[int(node.GetValue(Kratos.DENSITY))])

    def test_MapConditionContainerVariableDataHolderToNodalVariableDataHolder(self):
        def __GenerateEntity(center: Kratos.Array3, model_part: Kratos.ModelPart):
            model_part.CreateNewCondition("SurfaceCondition3D4N", model_part.NumberOfConditions() + 1, [
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([-0.5, -0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([0.5, -0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([0.5, 0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([-0.5, 0.5, 0.0]), model_part)
            ], model_part.GetProperties()[1])

        model = Kratos.Model()
        model_part = model.CreateModelPart("test1")
        model_part.CreateNewProperties(1)
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        for i in range(9):
            __GenerateEntity(Kratos.Array3([i // 3, i % 3, 0]), model_part)

        for condition in model_part.Conditions:
            condition.SetValue(Kratos.VELOCITY, Kratos.Array3([condition.Id, condition.Id + 1, condition.Id + 3]))
            condition.SetValue(Kratos.PRESSURE, condition.Id + 4)

        condition_container = KratosOA.ConditionContainerVariableDataHolder(model_part)
        neighbour_conditions = KratosOA.NodalContainerVariableDataHolder(model_part)
        mapped_values = KratosOA.NodalContainerVariableDataHolder(model_part)

        KratosOA.ContainerVariableDataHolderUtils.ComputeNumberOfNeighbourEntities(neighbour_conditions, condition_container)
        neighbour_conditions.AssignDataToContainerVariable(Kratos.YOUNG_MODULUS)

        condition_container.ReadDataFromContainerVariable(Kratos.VELOCITY)
        KratosOA.ContainerVariableDataHolderUtils.MapContainerVariableDataHolderToNodalVariableDataHolder(mapped_values, condition_container, neighbour_conditions)
        mapped_values.AssignDataToContainerVariable(Kratos.VELOCITY)

        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.ACCELERATION, model_part.Nodes)
        for condition in model_part.Conditions:
            for node in condition.GetGeometry():
                node[Kratos.ACCELERATION] += condition.GetValue(Kratos.VELOCITY) / node.GetValue(Kratos.YOUNG_MODULUS)

        for node in model_part.Nodes:
            self.assertVectorAlmostEqual(node[Kratos.ACCELERATION], node[Kratos.VELOCITY])

        condition_container.ReadDataFromContainerVariable(Kratos.PRESSURE)
        KratosOA.ContainerVariableDataHolderUtils.MapContainerVariableDataHolderToNodalVariableDataHolder(mapped_values, condition_container, neighbour_conditions)
        mapped_values.AssignDataToContainerVariable(Kratos.PRESSURE)

        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.DENSITY, model_part.Nodes)
        for condition in model_part.Conditions:
            for node in condition.GetGeometry():
                node[Kratos.DENSITY] += condition.GetValue(Kratos.PRESSURE) / node.GetValue(Kratos.YOUNG_MODULUS)

        for node in model_part.Nodes:
            self.assertEqual(node[Kratos.DENSITY], node[Kratos.PRESSURE])

    def test_MapElementContainerVariableDataHolderToNodalVariableDataHolder(self):
        def __GenerateEntity(center: Kratos.Array3, model_part: Kratos.ModelPart):
            model_part.CreateNewElement("Element2D4N", model_part.NumberOfElements() + 1, [
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([-0.5, -0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([0.5, -0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([0.5, 0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([-0.5, 0.5, 0.0]), model_part)
            ], model_part.GetProperties()[1])

        model = Kratos.Model()
        model_part = model.CreateModelPart("test1")
        model_part.CreateNewProperties(1)
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        for i in range(9):
            __GenerateEntity(Kratos.Array3([i // 3, i % 3, 0]), model_part)

        for element in model_part.Elements:
            element.SetValue(Kratos.VELOCITY, Kratos.Array3([element.Id, element.Id + 1, element.Id + 3]))
            element.SetValue(Kratos.PRESSURE, element.Id + 4)

        element_container = KratosOA.ElementContainerVariableDataHolder(model_part)
        neighbour_conditions = KratosOA.NodalContainerVariableDataHolder(model_part)
        mapped_values = KratosOA.NodalContainerVariableDataHolder(model_part)

        KratosOA.ContainerVariableDataHolderUtils.ComputeNumberOfNeighbourEntities(neighbour_conditions, element_container)
        neighbour_conditions.AssignDataToContainerVariable(Kratos.YOUNG_MODULUS)

        element_container.ReadDataFromContainerVariable(Kratos.VELOCITY)
        KratosOA.ContainerVariableDataHolderUtils.MapContainerVariableDataHolderToNodalVariableDataHolder(mapped_values, element_container, neighbour_conditions)
        mapped_values.AssignDataToContainerVariable(Kratos.VELOCITY)

        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.ACCELERATION, model_part.Nodes)
        for element in model_part.Elements:
            for node in element.GetGeometry():
                node[Kratos.ACCELERATION] += element.GetValue(Kratos.VELOCITY) / node.GetValue(Kratos.YOUNG_MODULUS)

        for node in model_part.Nodes:
            self.assertVectorAlmostEqual(node[Kratos.ACCELERATION], node[Kratos.VELOCITY])

        element_container.ReadDataFromContainerVariable(Kratos.PRESSURE)
        KratosOA.ContainerVariableDataHolderUtils.MapContainerVariableDataHolderToNodalVariableDataHolder(mapped_values, element_container, neighbour_conditions)
        mapped_values.AssignDataToContainerVariable(Kratos.PRESSURE)

        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.DENSITY, model_part.Nodes)
        for element in model_part.Elements:
            for node in element.GetGeometry():
                node[Kratos.DENSITY] += element.GetValue(Kratos.PRESSURE) / node.GetValue(Kratos.YOUNG_MODULUS)

        for node in model_part.Nodes:
            self.assertEqual(node[Kratos.DENSITY], node[Kratos.PRESSURE])

    def test_MapNodalVariableDataHolderToConditionContainerVariableDataHolder(self):
        def __GenerateEntity(center: Kratos.Array3, model_part: Kratos.ModelPart):
            model_part.CreateNewCondition("SurfaceCondition3D4N", model_part.NumberOfConditions() + 1, [
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([-0.5, -0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([0.5, -0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([0.5, 0.5, 0.0]), model_part),
                    TestContainerVariableDataHolderUtils.__GetNodeIndex(center + Kratos.Array3([-0.5, 0.5, 0.0]), model_part)
            ], model_part.GetProperties()[1])

        model = Kratos.Model()
        model_part = model.CreateModelPart("test1")
        model_part.CreateNewProperties(1)
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        for i in range(9):
            __GenerateEntity(Kratos.Array3([i // 3, i % 3, 0]), model_part)

        for node in model_part.Nodes:
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([node.Id, node.Id + 1, node.Id + 3]))
            node.SetValue(Kratos.PRESSURE, node.Id + 4)

        nodal_container = KratosOA.NodalContainerVariableDataHolder(model_part)
        mapped_value = KratosOA.ConditionContainerVariableDataHolder(model_part)

        nodal_container.ReadDataFromContainerVariable(Kratos.VELOCITY)
        KratosOA.ContainerVariableDataHolderUtils.MapNodalVariableDataHolderToContainerVariableDataHolder(mapped_value, nodal_container)
        mapped_value.AssignDataToContainerVariable(Kratos.ACCELERATION)

        for condition in model_part.Conditions:
            v = Kratos.Array3([0, 0, 0])
            for node in condition.GetGeometry():
                v += node.GetValue(Kratos.VELOCITY)
            self.assertVectorAlmostEqual(v / 4.0, condition.GetValue(Kratos.ACCELERATION))

        nodal_container.ReadDataFromContainerVariable(Kratos.PRESSURE)
        KratosOA.ContainerVariableDataHolderUtils.MapNodalVariableDataHolderToContainerVariableDataHolder(mapped_value, nodal_container)
        mapped_value.AssignDataToContainerVariable(Kratos.DENSITY)

        for condition in model_part.Conditions:
            v = 0.0
            for node in condition.GetGeometry():
                v += node.GetValue(Kratos.PRESSURE)
            self.assertEqual(v / 4.0, condition.GetValue(Kratos.DENSITY))

    def test_ComputeVariableDataHolderProductWithEntityMatrix(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("test1")
        model_part.AddNodalSolutionStepVariable(KratosOA.HELMHOLTZ_VAR_DENSITY)

        model_part.CreateNewNode(1, 0, 0, 0)
        model_part.CreateNewNode(2, 0, 1, 0)
        model_part.CreateNewNode(3, 1, 0, 0)
        model_part.CreateNewNode(4, 0, 0, 1)
        model_part.CreateNewNode(5, 1, 1, 0.5)

        properties = model_part.CreateNewProperties(1)
        model_part.CreateNewElement("HelmholtzBulkTopology3D4N", 1, [1, 2, 3, 4], properties)
        model_part.CreateNewElement("HelmholtzBulkTopology3D4N", 2, [2, 3, 4,5], properties)
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        for node in model_part.Nodes:
            node.SetValue(Kratos.PRESSURE, node.Id + 1)
            node.SetSolutionStepValue(KratosOA.HELMHOLTZ_VAR_DENSITY, node.Id + 1)

        nodal_values = KratosOA.NodalContainerVariableDataHolder(model_part)
        nodal_values.ReadDataFromContainerVariable(Kratos.PRESSURE)

        output_values = KratosOA.NodalContainerVariableDataHolder(model_part)
        KratosOA.ContainerVariableDataHolderUtils.ComputeVariableDataHolderProductWithEntityMatrix(output_values, nodal_values, KratosOA.HELMHOLTZ_MASS_MATRIX, model_part.Elements)

        # analytical calculation
        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.DENSITY, model_part.Nodes)
        element: Kratos.Element
        for element in model_part.Elements:
            v_in = Kratos.Vector(4)
            for i, node in enumerate(element.GetGeometry()):
                v_in[i] = node.GetSolutionStepValue(KratosOA.HELMHOLTZ_VAR_DENSITY)

            m = element.Calculate(KratosOA.HELMHOLTZ_MASS_MATRIX, model_part.ProcessInfo)
            v_out = m * v_in

            for i, node in enumerate(element.GetGeometry()):
                node.SetValue(Kratos.DENSITY, node.GetValue(Kratos.DENSITY) + v_out[i])

        analytical_values = KratosOA.NodalContainerVariableDataHolder(model_part)
        analytical_values.ReadDataFromContainerVariable(Kratos.DENSITY)

        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormL2(analytical_values - output_values), 0.0)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()