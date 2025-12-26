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
        cls.model_part.AddNodalSolutionStepVariable(Kratos.THICKNESS)
        with kratos_unittest.WorkFolderScope(".", __file__, True):
            ReadModelPart("model_part_utils_test/extruded_tetrahedra", cls.model_part)

        for node in cls.model_part.Nodes:
            id = node.Id
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetSolutionStepValue(Kratos.PRESSURE, id+3)
            node.SetSolutionStepValue(Kratos.DENSITY, id+4)
            node.SetSolutionStepValue(Kratos.THICKNESS, -id)

    # def test_ContainerVariableDataEntityMaxNormL2(self):
    #     a = Kratos.Expression.NodalExpression(self.model_part)

    #     Kratos.Expression.VariableExpressionIO.Read(a, Kratos.PRESSURE, True)
    #     self.assertEqual(KratosOA.ExpressionUtils.EntityMaxNormL2(a), 28)

    #     Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, True)
    #     self.assertEqual(KratosOA.ExpressionUtils.EntityMaxNormL2(a), math.sqrt(28**2 + 29**2 + 30**2))

    # def test_ContainerVariableDataInnerProduct(self):
    #     a = Kratos.Expression.NodalExpression(self.model_part)
    #     b = Kratos.Expression.NodalExpression(self.model_part)

    #     Kratos.Expression.VariableExpressionIO.Read(a, Kratos.PRESSURE, True)
    #     Kratos.Expression.VariableExpressionIO.Read(b, Kratos.DENSITY, True)

    #     self.assertEqual(Kratos.Expression.Utils.InnerProduct(a, b), 8100.0)

    # def test_ComputeNumberOfNeighbourConditions(self):
    #     neighbour_conditions = Kratos.Expression.NodalExpression(self.model_part)
    #     KratosOA.ExpressionUtils.ComputeNumberOfNeighbourConditions(neighbour_conditions)
    #     Kratos.Expression.VariableExpressionIO.Write(neighbour_conditions, Kratos.DENSITY, False)

    #     neighbour_map = {
    #         2.0: [1, 2, 4, 5, 7, 9, 10, 12, 15, 17, 20, 21, 22, 23, 24, 25],
    #         0.0: [3, 6, 8, 11, 13, 14, 16, 18, 19]
    #     }

    #     for node in self.model_part.Nodes:
    #         self.assertTrue(node.Id in neighbour_map[int(node.GetValue(Kratos.DENSITY))])

    # def test_ComputeNumberOfNeighbourElements(self):
    #     neighbour_elements = Kratos.Expression.NodalExpression(self.model_part)
    #     KratosOA.ExpressionUtils.ComputeNumberOfNeighbourElements(neighbour_elements)
    #     Kratos.Expression.VariableExpressionIO.Write(neighbour_elements, Kratos.DENSITY, False)

    #     neighbour_map = {
    #         1.0: [1, 9, 22, 25],
    #         2.0: [2, 4, 5, 7, 10, 12, 15, 17, 20, 21, 23, 24],
    #         4.0: [3, 6, 8, 11, 13, 14, 16, 18, 19]
    #     }

    #     for node in self.model_part.Nodes:
    #         self.assertTrue(node.Id in neighbour_map[int(node.GetValue(Kratos.DENSITY))])

    # def test_MapConditionVariableToNodalVariable(self):
    #     communicator: Kratos.Communicator = self.model_part.GetCommunicator()

    #     for condition in self.model_part.Conditions:
    #         condition.SetValue(Kratos.VELOCITY, Kratos.Array3([condition.Id, condition.Id + 1, condition.Id + 3]))
    #         condition.SetValue(Kratos.PRESSURE, condition.Id + 4)

    #     condition_container = Kratos.Expression.ConditionExpression(self.model_part)
    #     neighbour_conditions = Kratos.Expression.NodalExpression(self.model_part)
    #     mapped_values = Kratos.Expression.NodalExpression(self.model_part)

    #     KratosOA.ExpressionUtils.ComputeNumberOfNeighbourConditions(neighbour_conditions)
    #     Kratos.Expression.VariableExpressionIO.Write(neighbour_conditions, Kratos.YOUNG_MODULUS, False)

    #     Kratos.Expression.VariableExpressionIO.Read(condition_container, Kratos.VELOCITY)
    #     KratosOA.ExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, condition_container, neighbour_conditions)
    #     Kratos.Expression.VariableExpressionIO.Write(mapped_values, Kratos.VELOCITY, False)

    #     Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.ACCELERATION, self.model_part.Nodes)
    #     for condition in self.model_part.Conditions:
    #         for node in condition.GetGeometry():
    #             node[Kratos.ACCELERATION] += condition.GetValue(Kratos.VELOCITY) / node.GetValue(Kratos.YOUNG_MODULUS)

    #     communicator.AssembleNonHistoricalData(Kratos.ACCELERATION)

    #     for node in self.model_part.Nodes:
    #         self.assertVectorAlmostEqual(node[Kratos.ACCELERATION], node[Kratos.VELOCITY])

    #     Kratos.Expression.VariableExpressionIO.Read(condition_container, Kratos.PRESSURE)
    #     KratosOA.ExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, condition_container, neighbour_conditions)
    #     Kratos.Expression.VariableExpressionIO.Write(mapped_values, Kratos.PRESSURE, False)

    #     Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.DENSITY, self.model_part.Nodes)
    #     for condition in self.model_part.Conditions:
    #         for node in condition.GetGeometry():
    #             node[Kratos.DENSITY] += condition.GetValue(Kratos.PRESSURE) / node.GetValue(Kratos.YOUNG_MODULUS)

    #     communicator.AssembleNonHistoricalData(Kratos.DENSITY)

    #     for node in self.model_part.Nodes:
    #         self.assertAlmostEqual(node[Kratos.DENSITY], node[Kratos.PRESSURE], 9)

    # def test_MapElementVariableToNodalVariable(self):
    #     communicator: Kratos.Communicator = self.model_part.GetCommunicator()

    #     for element in self.model_part.Elements:
    #         element.SetValue(Kratos.VELOCITY, Kratos.Array3([element.Id, element.Id + 1, element.Id + 3]))
    #         element.SetValue(Kratos.PRESSURE, element.Id + 4)

    #     element_container = Kratos.Expression.ElementExpression(self.model_part)
    #     neighbour_elements = Kratos.Expression.NodalExpression(self.model_part)
    #     mapped_values = Kratos.Expression.NodalExpression(self.model_part)

    #     KratosOA.ExpressionUtils.ComputeNumberOfNeighbourElements(neighbour_elements)
    #     Kratos.Expression.VariableExpressionIO.Write(neighbour_elements, Kratos.YOUNG_MODULUS, False)

    #     Kratos.Expression.VariableExpressionIO.Read(element_container, Kratos.VELOCITY)
    #     KratosOA.ExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, element_container, neighbour_elements)
    #     Kratos.Expression.VariableExpressionIO.Write(mapped_values, Kratos.VELOCITY, False)

    #     Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.ACCELERATION, self.model_part.Nodes)
    #     for element in self.model_part.Elements:
    #         for node in element.GetGeometry():
    #             node[Kratos.ACCELERATION] += element.GetValue(Kratos.VELOCITY) / node.GetValue(Kratos.YOUNG_MODULUS)

    #     communicator.AssembleNonHistoricalData(Kratos.ACCELERATION)

    #     for node in self.model_part.Nodes:
    #         self.assertVectorAlmostEqual(node[Kratos.ACCELERATION], node[Kratos.VELOCITY])

    #     Kratos.Expression.VariableExpressionIO.Read(element_container, Kratos.PRESSURE)
    #     KratosOA.ExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, element_container, neighbour_elements)
    #     Kratos.Expression.VariableExpressionIO.Write(mapped_values, Kratos.PRESSURE, False)

    #     Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.DENSITY, self.model_part.Nodes)
    #     for element in self.model_part.Elements:
    #         for node in element.GetGeometry():
    #             node[Kratos.DENSITY] += element.GetValue(Kratos.PRESSURE) / node.GetValue(Kratos.YOUNG_MODULUS)

    #     communicator.AssembleNonHistoricalData(Kratos.DENSITY)

    #     for node in self.model_part.Nodes:
    #         self.assertEqual(node[Kratos.DENSITY], node[Kratos.PRESSURE])

    # def test_MapNodalVariableToConditionVariable(self):
    #     for node in self.model_part.Nodes:
    #         node.SetValue(Kratos.VELOCITY, Kratos.Array3([node.Id, node.Id + 1, node.Id + 3]))
    #         node.SetValue(Kratos.PRESSURE, node.Id + 4)

    #     nodal_container = Kratos.Expression.NodalExpression(self.model_part)
    #     mapped_value = Kratos.Expression.ConditionExpression(self.model_part)

    #     Kratos.Expression.VariableExpressionIO.Read(nodal_container, Kratos.VELOCITY, False)
    #     KratosOA.ExpressionUtils.MapNodalVariableToContainerVariable(mapped_value, nodal_container)
    #     Kratos.Expression.VariableExpressionIO.Write(mapped_value, Kratos.ACCELERATION)

    #     for condition in self.model_part.Conditions:
    #         v = Kratos.Array3([0, 0, 0])
    #         for node in condition.GetGeometry():
    #             v += node.GetValue(Kratos.VELOCITY)
    #         self.assertVectorAlmostEqual(v / 2.0, condition.GetValue(Kratos.ACCELERATION))

    #     Kratos.Expression.VariableExpressionIO.Read(nodal_container, Kratos.PRESSURE, False)
    #     KratosOA.ExpressionUtils.MapNodalVariableToContainerVariable(mapped_value, nodal_container)
    #     Kratos.Expression.VariableExpressionIO.Write(mapped_value, Kratos.DENSITY)

    #     for condition in self.model_part.Conditions:
    #         v = 0.0
    #         for node in condition.GetGeometry():
    #             v += node.GetValue(Kratos.PRESSURE)
    #         self.assertEqual(v / 2.0, condition.GetValue(Kratos.DENSITY))

    # def test_ProductWithEntityMatrix(self):
    #     if (IsDistributedRun()):
    #         self.skipTest("Skipping since ProductWithEntityMatrix does not support MPI yet.")

    #     number_of_nodes = self.model_part.NumberOfNodes()

    #     a = Kratos.Expression.NodalExpression(self.model_part)
    #     Kratos.Expression.VariableExpressionIO.Read(a, Kratos.PRESSURE, True)

    #     m = Kratos.Matrix(number_of_nodes, number_of_nodes)
    #     for i in range(number_of_nodes):
    #         for j in range(number_of_nodes):
    #             m[i, j] = (i + 1) * (j + 1)

    #     b = Kratos.Expression.NodalExpression(self.model_part)
    #     KratosOA.ExpressionUtils.ProductWithEntityMatrix(b, m, a)
    #     Kratos.Expression.VariableExpressionIO.Write(b, Kratos.DENSITY, True)

    #     for i, node_b in enumerate(b.GetContainer()):
    #         v = 0
    #         for j, node_a in enumerate(a.GetContainer()):
    #             v += m[i, j] * node_a.GetSolutionStepValue(Kratos.PRESSURE)
    #         self.assertEqual(v, node_b.GetSolutionStepValue(Kratos.DENSITY))

    # def test_ProductWithEntityMatrixSparse(self):
    #     if (IsDistributedRun()):
    #         self.skipTest("Skipping since ProductWithEntityMatrix does not support MPI yet.")

    #     number_of_nodes = self.model_part.NumberOfNodes()

    #     a = Kratos.Expression.NodalExpression(self.model_part)
    #     Kratos.Expression.VariableExpressionIO.Read(a, Kratos.PRESSURE, True)

    #     dense_m = Kratos.Matrix(number_of_nodes, number_of_nodes)
    #     for i in range(number_of_nodes):
    #         for j in range(number_of_nodes):
    #             dense_m[i, j] = int((i + 1) * (j + 1) % 10)

    #     # now add values to the sparse matrix
    #     sparse_m = Kratos.CompressedMatrix(number_of_nodes, number_of_nodes)
    #     for i in range(number_of_nodes):
    #         for j in range(number_of_nodes):
    #             if dense_m[i, j] != 0.0:
    #                 sparse_m[i, j] = dense_m[i, j]

    #     dense_b = Kratos.Expression.NodalExpression(self.model_part)
    #     KratosOA.ExpressionUtils.ProductWithEntityMatrix(dense_b, dense_m, a)

    #     sparse_b = Kratos.Expression.NodalExpression(self.model_part)
    #     KratosOA.ExpressionUtils.ProductWithEntityMatrix(sparse_b, sparse_m, a)

    #     self.assertEqual(Kratos.Expression.Utils.InnerProduct(dense_b - sparse_b, dense_b - sparse_b), 0)

    # def test_Transpose(self):
    #     matrix_size = 10

    #     # first build the normal matrix
    #     dense_m = Kratos.Matrix(
    #         [
    #             [5, 0, 0, 2, 0, 0, 0, 0, 0],
    #             [0, 0, 3, 2, 3, 0, 0, 0, 0],
    #             [0, 2, 0, 2, 0, 2, 0, 6, 0],
    #             [0, 0, 0, 2, 0, 4, 0, 5, 0],
    #             [0, 0, 2, 0, 0, 2, 0, 0, 0],
    #             [0, 0, 0, 0, 0, 0, 0, 3, 3],
    #             [0, 2, 0, 0, 3, 0, 4, 0, 0],
    #             [3, 0, 0, 2, 0, 0, 0, 2, 0],
    #             [9, 0, 0, 0, 0, 7, 0, 7, 0],
    #             [0, 0, 7, 2, 0, 0, 0, 0, 0]
    #         ]
    #     )
    #     transpose_dense_m = Kratos.Matrix()
    #     KratosOA.ExpressionUtils.Transpose(transpose_dense_m, dense_m)
    #     for i in range(dense_m.Size1()):
    #         for j in range(dense_m.Size2()):
    #             self.assertEqual(transpose_dense_m[j, i], dense_m[i, j])

    #     # now add values to the sparse matrix
    #     sparse_m = Kratos.CompressedMatrix(matrix_size, matrix_size)
    #     for i in range(dense_m.Size1()):
    #         for j in range(dense_m.Size2()):
    #             if dense_m[i, j] != 0.0:
    #                 sparse_m[i, j] = dense_m[i, j]

    #     transpose_sparse_m = Kratos.CompressedMatrix()
    #     KratosOA.ExpressionUtils.Transpose(transpose_sparse_m, sparse_m)
    #     for i in range(dense_m.Size1()):
    #         for j in range(dense_m.Size2()):
    #             self.assertEqual(transpose_sparse_m[j, i], dense_m[i, j])

    # def test_ComputeVariableDataHolderProductWithEntityMatrix(self):
    #     for node in self.model_part.Nodes:
    #         node.SetValue(Kratos.PRESSURE, node.Id + 1)
    #         node.SetSolutionStepValue(KratosOA.HELMHOLTZ_VAR_DENSITY, node.Id + 1)

    #     nodal_values = Kratos.Expression.NodalExpression(self.model_part)
    #     Kratos.Expression.VariableExpressionIO.Read(nodal_values, Kratos.PRESSURE, False)

    #     output_values = Kratos.Expression.NodalExpression(self.model_part)
    #     KratosOA.ExpressionUtils.ComputeNodalVariableProductWithEntityMatrix(output_values, nodal_values, KratosOA.HELMHOLTZ_MASS_MATRIX, self.model_part.Elements)

    #     # analytical calculation
    #     Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.DENSITY, self.model_part.Nodes)
    #     element: Kratos.Element
    #     for element in self.model_part.Elements:
    #         v_in = Kratos.Vector(4)
    #         for i, node in enumerate(element.GetGeometry()):
    #             v_in[i] = node.GetSolutionStepValue(KratosOA.HELMHOLTZ_VAR_DENSITY)

    #         m = element.Calculate(KratosOA.HELMHOLTZ_MASS_MATRIX, self.model_part.ProcessInfo)
    #         v_out = m * v_in

    #         for i, node in enumerate(element.GetGeometry()):
    #             node.SetValue(Kratos.DENSITY, node.GetValue(Kratos.DENSITY) + v_out[i])

    #     self.model_part.GetCommunicator().AssembleNonHistoricalData(Kratos.DENSITY)

    #     analytical_values = Kratos.Expression.NodalExpression(self.model_part)
    #     Kratos.Expression.VariableExpressionIO.Read(analytical_values, Kratos.DENSITY, False)

    #     self.assertEqual(Kratos.Expression.Utils.NormL2(analytical_values - output_values), 0.0)

    # def test_Scopes(self):
    #     a = Kratos.Expression.NodalExpression(self.model_part)
    #     b = Kratos.Expression.NodalExpression(self.model_part)

    #     Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, True)
    #     Kratos.Expression.VariableExpressionIO.Read(b, Kratos.PRESSURE, True)

    #     def func(a):
    #         c = a * 2
    #         d = c * 3
    #         return d

    #     c = func(a)
    #     d = func(b)

    #     e = c.Clone()
    #     Kratos.Expression.VariableExpressionIO.Write(e, Kratos.ACCELERATION, False)
    #     f = d.Clone()
    #     Kratos.Expression.VariableExpressionIO.Write(f, Kratos.DENSITY, False)

    #     for node in self.model_part.Nodes:
    #         self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.VELOCITY) * 6, node.GetValue(Kratos.ACCELERATION))
    #         self.assertEqual(node.GetSolutionStepValue(Kratos.PRESSURE) * 6, node.GetValue(Kratos.DENSITY))

    # def test_CollectiveExpressionsNormInf(self):
    #     a = Kratos.Expression.NodalExpression(self.model_part)
    #     b = Kratos.Expression.ElementExpression(self.model_part)

    #     Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, True)
    #     KratosOA.PropertiesVariableExpressionIO.Read(b, Kratos.PRESSURE)

    #     collective_1 = KratosOA.CollectiveExpression([a, b])
    #     self.assertEqual(KratosOA.ExpressionUtils.NormInf(collective_1), max(Kratos.Expression.Utils.NormInf(a), Kratos.Expression.Utils.NormInf(b)))

    # def test_CollectiveExpressionsNormL2(self):
    #     a = Kratos.Expression.NodalExpression(self.model_part)
    #     b = Kratos.Expression.ElementExpression(self.model_part)

    #     Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, True)
    #     KratosOA.PropertiesVariableExpressionIO.Read(b, Kratos.PRESSURE)

    #     collective_1 = KratosOA.CollectiveExpression([a, b])
    #     self.assertEqual(KratosOA.ExpressionUtils.NormL2(collective_1), math.sqrt(Kratos.Expression.Utils.NormL2(a)**2 + Kratos.Expression.Utils.NormL2(b)**2))

    # def test_CollectiveExpressionsInnerProduct(self):
    #     a = Kratos.Expression.NodalExpression(self.model_part)
    #     b = Kratos.Expression.ElementExpression(self.model_part)

    #     Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, True)
    #     KratosOA.PropertiesVariableExpressionIO.Read(b, Kratos.PRESSURE)

    #     collective_1 = KratosOA.CollectiveExpression([a, b])
    #     self.assertEqual(KratosOA.ExpressionUtils.InnerProduct(collective_1, collective_1), Kratos.Expression.Utils.InnerProduct(a, a) + Kratos.Expression.Utils.InnerProduct(b, b))

    def test_GetGradientExpression(self):
        input_nodal = Kratos.Expression.NodalExpression(self.model_part)
        output_elemental = Kratos.Expression.ElementExpression(self.model_part)
        dimension = 3

        Kratos.Expression.VariableExpressionIO.Read(input_nodal, Kratos.PRESSURE, False)

        KratosOA.ExpressionUtils.GetGradientExpression(output_elemental, input_nodal, dimension)
        
        # print(f"\nNODAL INPUT (GRAD): {input_nodal.Evaluate()}\nELEMENTAL OUTPUT (GRAD): {output_elemental.Evaluate()}\n")

    def test_ProjectElementalToNodalViaShapeFunctions(self):
        input_elemental = Kratos.Expression.ElementExpression(self.model_part)
        output_nodal = Kratos.Expression.NodalExpression(self.model_part)

        element: Kratos.Element
        for element in self.model_part.Elements:
            element.SetValue(Kratos.PRESSURE, element.Id + 4)

        Kratos.Expression.VariableExpressionIO.Read(input_elemental, Kratos.PRESSURE)

        KratosOA.ExpressionUtils.ProjectElementalToNodalViaShapeFunctions(output_nodal, input_elemental)
        
        # print(f"\nELEMENTAL INPUT: {input_elemental.Evaluate()}\nNODAL OUTPUT: {output_nodal.Evaluate()}\n")

    def test_ProjectNodalToElementalViaShapeFunctions(self):
        input_nodal = Kratos.Expression.NodalExpression(self.model_part)
        output_elemental = Kratos.Expression.ElementExpression(self.model_part)

        node: Kratos.Node
        for node in self.model_part.Nodes:
            node.SetValue(Kratos.DENSITY, node.Id + 4)

        Kratos.Expression.VariableExpressionIO.Read(input_nodal, Kratos.DENSITY, False)

        KratosOA.ExpressionUtils.ProjectNodalToElementalViaShapeFunctions(output_elemental, input_nodal)
        
        # print(f"\nNODAL INPUT: {input_nodal.Evaluate()}\nELEMENTAL OUTPUT: {output_elemental.Evaluate()}\n")

    def test_HamilotinanUpdate(self):
        input_phi = Kratos.Expression.NodalExpression(self.model_part)
        input_v = Kratos.Expression.ElementExpression(self.model_part)
        input_grad = Kratos.Expression.ElementExpression(self.model_part)
        output_phi = Kratos.Expression.NodalExpression(self.model_part)

        node: Kratos.Node
        for node in self.model_part.Nodes:
            node.SetValue(Kratos.PRESSURE, node.Id + 4)

        Kratos.Expression.VariableExpressionIO.Read(input_phi, Kratos.PRESSURE, False)
        element: Kratos.Element
        for element in self.model_part.Elements:
            element.SetValue(Kratos.DENSITY, element.Id + element.Id*0.2)
            element.SetValue(Kratos.YOUNG_MODULUS, element.Id*0.1)
            
        Kratos.Expression.VariableExpressionIO.Read(input_v, Kratos.DENSITY)
        Kratos.Expression.VariableExpressionIO.Read(input_grad, Kratos.YOUNG_MODULUS)

        KratosOA.ExpressionUtils.HamilotinanUpdate(input_phi, output_phi, input_v, input_grad)
        # print(f"INPUT_PHI: {input_phi.Evaluate()}\nINPUT_V: {input_v.Evaluate()}\nINPUT_GRAD: {input_grad.Evaluate()}\nOUTPUT PHI: {output_phi.Evaluate()}")

    def test_Heaviside(self):
        import csv
        import math
        from pathlib import Path

        # Heaviside parameters to sweep
        precisions = [50, 100, 200, 1000, 2000]
        k = 10.0

        # Expression containers
        input_elemental = Kratos.Expression.ElementExpression(self.model_part)

        # --- Build synthetic phi in [-1, 1] over the elements ---
        elements = sorted(self.model_part.Elements, key=lambda e: e.Id)
        n_elems = len(elements)

        element_ids = []
        phi_vals = []

        for i, element in enumerate(elements):
            # linearly spaced in [-1, 1]
            if n_elems == 1:
                phi = 0.0
            else:
                phi = -1.0 + 2.0 * i / (n_elems - 1)

            element.SetValue(Kratos.DENSITY, phi)
            element_ids.append(element.Id)
            phi_vals.append(phi)

        # Read DENSITY into the input expression
        Kratos.Expression.VariableExpressionIO.Read(input_elemental, Kratos.DENSITY)

        # --- Analytic reference H_ref(φ) computed ONCE ---
        ref_vals = [
            1.0 / (1.0 + math.exp(-2.0 * k * phi))
            for phi in phi_vals
        ]

        # --- Evaluate Kratos Heaviside for each precision ---
        kratos_results = {}  # precision -> list of H values

        for precision in precisions:
            output_elemental = Kratos.Expression.ElementExpression(self.model_part)
            KratosOA.ExpressionUtils.Heaviside(output_elemental, input_elemental, precision, k)
            kratos_results[precision] = list(output_elemental.Evaluate())

        # --- Optional check: compare one precision against reference ---
        # Use a reasonable decimal accuracy; 14–15 is usually plenty for doubles.
        check_precision = precisions[-1]  # e.g. 2000
        for i in range(n_elems):
            self.assertAlmostEqual(
                kratos_results[check_precision][i],
                ref_vals[i],
                15  # you can bump this to 17 if you want to reproduce your experiment
            )

        # --- Write everything to CSV once ---
        out_dir = Path("analysis")
        out_dir.mkdir(exist_ok=True)
        out_path =  Path("heaviside_test.csv")

        with out_path.open("w", newline="") as f:
            writer = csv.writer(f)
            header = ["element_id", "phi", "H_ref"] + [f"H_prec_{p}" for p in precisions]
            writer.writerow(header)

            for idx, eid in enumerate(element_ids):
                row = [eid, phi_vals[idx], ref_vals[idx]]
                for p in precisions:
                    row.append(kratos_results[p][idx])
                writer.writerow(row)

        print(f"Heaviside test data written to {out_path}")

    def test_DiracDelta(self):
        import csv
        import math
        from pathlib import Path

        # Precisions to test (same style as for Heaviside)
        precisions = [50, 100, 200, 1000]
        k = 10.0

        # ContainerExpressions
        input_elemental = Kratos.Expression.ElementExpression(self.model_part)

        # --- Build synthetic phi in [-1, 1] over the elements ---
        elements = sorted(self.model_part.Elements, key=lambda e: e.Id)
        n_elems = len(elements)

        element_ids = []
        phi_vals = []

        for i, element in enumerate(elements):
            if n_elems == 1:
                phi = 0.0
            else:
                phi = -1.0 + 2.0 * i / (n_elems - 1)  # linear in [-1,1]

            element.SetValue(Kratos.DENSITY, phi)
            element_ids.append(element.Id)
            phi_vals.append(phi)

        # Read DENSITY into the input expression
        Kratos.Expression.VariableExpressionIO.Read(input_elemental, Kratos.DENSITY)

        # --- Analytic reference δ_ref(φ) = dH/dφ ---
        delta_ref = []
        for phi in phi_vals:
            s = -2.0 * k * phi
            e = math.exp(s)
            # H' = 2k * e^{-2k phi} / (1 + e^{-2k phi})^2
            delta_ref.append(2.0 * k * e / (1.0 + e) ** 2)

        # --- Evaluate Kratos DiracDelta for each precision ---
        kratos_results = {}  # precision -> list of δ values

        for precision in precisions:
            output_elemental = Kratos.Expression.ElementExpression(self.model_part)
            KratosOA.ExpressionUtils.DiracDelta(
                output_elemental, input_elemental, precision, k
            )
            kratos_results[precision] = list(output_elemental.Evaluate())

        # --- Optional consistency check against analytic reference ---
        # Use a realistic decimal accuracy for doubles (12–15).
        check_precision = precisions[-1]  # e.g. 1000
        for i in range(n_elems):
            self.assertAlmostEqual(
                kratos_results[check_precision][i],
                delta_ref[i],
                12  # tighten/loosen as you like
            )

        # --- Write everything to CSV ---
        out_dir = Path("analysis")
        out_dir.mkdir(exist_ok=True)
        out_path = out_dir / "dirac_test.csv"

        with out_path.open("w", newline="") as f:
            writer = csv.writer(f)
            header = ["element_id", "phi", "delta_ref"] + [f"delta_prec_{p}" for p in precisions]
            writer.writerow(header)

            for idx, eid in enumerate(element_ids):
                row = [eid, phi_vals[idx], delta_ref[idx]]
                for p in precisions:
                    row.append(kratos_results[p][idx])
                writer.writerow(row)

        print(f"Dirac delta test data written to {out_path}")




if __name__ == "__main__":
    kratos_unittest.main()