
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestExtractBoundaryNodesProcess(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()

    @staticmethod
    def __GetNodeIndex(model_part: Kratos.ModelPart, position: Kratos.Array3):
        for node in model_part.Nodes:
            distance = (position[0] - node.X)**2 + (position[1] - node.Y)**2 + (position[2] - node.Z)**2
            if distance < 1e-8:
                return node.Id

        node = model_part.CreateNewNode(model_part.NumberOfNodes(), position[0], position[1], position[2])
        return node.Id

    def test_ExtractBoundaryNodesProcess2D(self):
        model_part = self.model.CreateModelPart("test_surfaces")
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
        properties = model_part.CreateNewProperties(0)

        grid_size = 3
        for i in range(9):
            center = Kratos.Array3([i // grid_size, i % grid_size, 0])
            model_part.CreateNewElement("Element2D4N", i+1, [
                TestExtractBoundaryNodesProcess.__GetNodeIndex(model_part, center + Kratos.Array3([-0.5, -0.5, 0.0])),
                TestExtractBoundaryNodesProcess.__GetNodeIndex(model_part, center + Kratos.Array3([0.5, -0.5, 0.0])),
                TestExtractBoundaryNodesProcess.__GetNodeIndex(model_part, center + Kratos.Array3([0.5, 0.5, 0.0])),
                TestExtractBoundaryNodesProcess.__GetNodeIndex(model_part, center + Kratos.Array3([-0.5, 0.5, 0.0]))
            ], properties)


        process = KratosOA.ExtractBoundaryNodesProcess(self.model, Kratos.Parameters("""{
            "model_part_name": "test_surfaces",
            "boundary_sub_model_part_name": "Auto_Boundary"
        }"""))
        process.Check()
        process.ExecuteInitialize()

        self.assertTrue(model_part.HasSubModelPart("Auto_Boundary"))

        boundary_model_part = model_part.GetSubModelPart("Auto_Boundary")
        self.assertEqual(boundary_model_part.NumberOfNodes(), 12)

        number_of_nodes_left = 0
        number_of_nodes_right = 0
        number_of_nodes_top = 0
        number_of_nodes_bottom = 0
        for node in boundary_model_part.Nodes:
            if node.X == -0.5:
                number_of_nodes_left += 1
            if node.X == 2.5:
                number_of_nodes_right += 1
            if node.Y == -0.5:
                number_of_nodes_top += 1
            if node.Y == 2.5:
                number_of_nodes_bottom += 1

        self.assertEqual(number_of_nodes_left, 4)
        self.assertEqual(number_of_nodes_right, 4)
        self.assertEqual(number_of_nodes_top, 4)
        self.assertEqual(number_of_nodes_bottom, 4)

    def test_ExtractBoundaryNodesProcess3D(self):
        model_part = self.model.CreateModelPart("test_solids")
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        properties = model_part.CreateNewProperties(0)

        grid_size = 3
        for i in range(27):
            center = Kratos.Array3([i % grid_size, (i // grid_size) % grid_size, i // (grid_size ** 2)])
            model_part.CreateNewElement("Element3D8N", i+1, [
                TestExtractBoundaryNodesProcess.__GetNodeIndex(model_part, center + Kratos.Array3([-0.5, -0.5, -0.5])),
                TestExtractBoundaryNodesProcess.__GetNodeIndex(model_part, center + Kratos.Array3([0.5, -0.5, -0.5])),
                TestExtractBoundaryNodesProcess.__GetNodeIndex(model_part, center + Kratos.Array3([0.5, 0.5, -0.5])),
                TestExtractBoundaryNodesProcess.__GetNodeIndex(model_part, center + Kratos.Array3([-0.5, 0.5, -0.5])),
                TestExtractBoundaryNodesProcess.__GetNodeIndex(model_part, center + Kratos.Array3([-0.5, -0.5, 0.5])),
                TestExtractBoundaryNodesProcess.__GetNodeIndex(model_part, center + Kratos.Array3([0.5, -0.5, 0.5])),
                TestExtractBoundaryNodesProcess.__GetNodeIndex(model_part, center + Kratos.Array3([0.5, 0.5, 0.5])),
                TestExtractBoundaryNodesProcess.__GetNodeIndex(model_part, center + Kratos.Array3([-0.5, 0.5, 0.5]))
            ], properties)

        process = KratosOA.ExtractBoundaryNodesProcess(self.model, Kratos.Parameters("""{
            "model_part_name": "test_solids",
            "boundary_sub_model_part_name": "Auto_Boundary"
        }"""))
        process.Check()
        process.ExecuteInitialize()

        self.assertTrue(model_part.HasSubModelPart("Auto_Boundary"))

        boundary_model_part = model_part.GetSubModelPart("Auto_Boundary")
        self.assertEqual(boundary_model_part.NumberOfNodes(), 16+2*12+8+4*2)

        number_of_nodes_left = 0
        number_of_nodes_right = 0
        number_of_nodes_top = 0
        number_of_nodes_bottom = 0
        number_of_side_1 = 0
        number_of_side_2 = 0
        for node in boundary_model_part.Nodes:
            if node.X == -0.5:
                number_of_nodes_left += 1
            if node.X == 2.5:
                number_of_nodes_right += 1
            if node.Y == -0.5:
                number_of_nodes_top += 1
            if node.Y == 2.5:
                number_of_nodes_bottom += 1
            if node.Z == -0.5:
                number_of_side_1 += 1
            if node.Z == 2.5:
                number_of_side_2 += 1

        self.assertEqual(number_of_nodes_left, 16)
        self.assertEqual(number_of_nodes_right, 16)
        self.assertEqual(number_of_nodes_top, 16)
        self.assertEqual(number_of_nodes_bottom, 16)
        self.assertEqual(number_of_side_1, 16)
        self.assertEqual(number_of_side_2, 16)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()