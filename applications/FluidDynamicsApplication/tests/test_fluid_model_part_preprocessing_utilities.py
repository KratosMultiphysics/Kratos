import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestFluidModelPartPreprocessingUtilities(KratosUnittest.TestCase):
    def test_CreateModelPartForCommenInterface(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("test")

        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
        model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        model_part.CreateNewNode(4, 1.0, 0.0, 0.0)

        wall_model_part = model_part.CreateSubModelPart("wall")
        wall_model_part.AddNodes([1, 2, 4])
        side_model_part = model_part.CreateSubModelPart("side")
        side_model_part.AddNodes([1, 2])
        none_model_part = model_part.CreateSubModelPart("none")
        none_model_part.AddNodes([3, 4])

        self.assertTrue(KratosFluid.FluidModelPartPreProcessingUtilities.CreateModelPartForCommenInterface(model_part, "side_wall", ["wall", "side"]))
        self.assertFalse(KratosFluid.FluidModelPartPreProcessingUtilities.CreateModelPartForCommenInterface(model_part, "side_none", ["none", "side"]))

        for check_node in model["test.side_wall"].Nodes:
            self.assertTrue((check_node.Id in [node.Id for node in wall_model_part.Nodes]) and (check_node.Id in [node.Id for node in side_model_part.Nodes]))

        self.assertEqual(model["test.side_none"].NumberOfNodes(), 0)

    def test_GetElementIdsWithAllNodesOnBoundaries(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2

        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
        model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        model_part.CreateNewNode(4, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(5, 2.0, 1.0, 0.0)
        model_part.CreateNewNode(6, 2.0, 0.0, 0.0)

        wall_model_part = model_part.CreateSubModelPart("wall")
        wall_model_part.AddNodes([1, 2, 4])

        side_model_part = model_part.CreateSubModelPart("side")
        side_model_part.AddNodes([4, 6, 5])

        properties = model_part.GetProperties()[0]
        model_part.CreateNewElement("Element2D3N", 1, [1, 2, 4], properties)
        model_part.CreateNewElement("Element2D3N", 2, [2, 3, 4], properties)
        model_part.CreateNewElement("Element2D3N", 3, [3, 5, 4], properties)
        model_part.CreateNewElement("Element2D3N", 4, [4, 5, 6], properties)

        element_ids = KratosFluid.FluidModelPartPreProcessingUtilities.GetElementIdsWithAllNodesOnBoundaries(model_part, ["wall", "side"])
        self.assertEqual(element_ids, [1, 4])

    def test_FixBoundaryElementsTriangle(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2

        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
        model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        model_part.CreateNewNode(4, 1.0, 0.0, 0.0)

        wall_model_part = model_part.CreateSubModelPart("wall")
        wall_model_part.AddNodes([1, 2, 4])

        properties = model_part.GetProperties()[0]
        model_part.CreateNewElement("Element2D3N", 1, [1, 2, 4], properties)
        model_part.CreateNewElement("Element2D3N", 2, [2, 3, 4], properties)

        element_ids = KratosFluid.FluidModelPartPreProcessingUtilities.GetElementIdsWithAllNodesOnBoundaries(model_part, ["wall"])
        self.assertEqual(element_ids, [1])
        KratosFluid.FluidModelPartPreProcessingUtilities.BreakElements(model_part, "Element2D3N", element_ids)
        element_ids = KratosFluid.FluidModelPartPreProcessingUtilities.GetElementIdsWithAllNodesOnBoundaries(model_part, ["wall"])
        self.assertEqual(element_ids, [])
        # Kratos.ModelPartIO("test", Kratos.IO.IGNORE_VARIABLES_ERROR | Kratos.IO.WRITE | Kratos.IO.MESH_ONLY).WriteModelPart(model_part)

    def test_FixBoundaryElementsTetrahedra(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("test")

        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
        model_part.CreateNewNode(3, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(4, 0.5, 0.5, 1.0)
        model_part.CreateNewNode(5, 1.0, 1.0, 1.0)

        wall_model_part = model_part.CreateSubModelPart("wall")
        wall_model_part.AddNodes([1, 2, 3, 4])

        properties = model_part.GetProperties()[0]
        model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], properties)
        model_part.CreateNewElement("Element3D4N", 2, [2, 3, 4, 5], properties)

        element_ids = KratosFluid.FluidModelPartPreProcessingUtilities.GetElementIdsWithAllNodesOnBoundaries(model_part, ["wall"])
        self.assertEqual(element_ids, [1])
        KratosFluid.FluidModelPartPreProcessingUtilities.BreakElements(model_part, "Element3D4N", element_ids)
        element_ids = KratosFluid.FluidModelPartPreProcessingUtilities.GetElementIdsWithAllNodesOnBoundaries(model_part, ["wall"])
        self.assertEqual(element_ids, [])
        # Kratos.ModelPartIO("test", Kratos.IO.IGNORE_VARIABLES_ERROR | Kratos.IO.WRITE | Kratos.IO.MESH_ONLY).WriteModelPart(model_part)

if __name__ == '__main__':
    KratosUnittest.main()