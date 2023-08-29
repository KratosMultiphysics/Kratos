import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.HDF5Application as HDF5Application


class TestVertex(KratosUnittest.TestCase):

    @staticmethod
    def MakeModelPart():
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("main")

        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        node = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [0.0, 0.0, 0.0])
        node.SetValue(KratosMultiphysics.DISPLACEMENT, [0.0, 0.0, 0.0])
        node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, 0.0])
        node.SetValue(KratosMultiphysics.REACTION, [0.0, 0.0, 0.0])

        node = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [1.0, 0.0, 0.0])
        node.SetValue(KratosMultiphysics.DISPLACEMENT, [1.0, 0.0, 0.0])
        node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, 1.0])
        node.SetValue(KratosMultiphysics.REACTION, [0.0, 0.0, 1.0])

        node = model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [2.0, 0.0, 0.0])
        node.SetValue(KratosMultiphysics.DISPLACEMENT, [2.0, 0.0, 0.0])
        node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, 2.0])
        node.SetValue(KratosMultiphysics.REACTION, [0.0, 0.0, 2.0])

        node = model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [3.0, 0.0, 0.0])
        node.SetValue(KratosMultiphysics.DISPLACEMENT, [3.0, 0.0, 0.0])
        node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, 3.0])
        node.SetValue(KratosMultiphysics.REACTION, [0.0, 0.0, 3.0])

        node = model_part.CreateNewNode(5, 2.0, 0.0, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [4.0, 0.0, 0.0])
        node.SetValue(KratosMultiphysics.DISPLACEMENT, [4.0, 0.0, 0.0])
        node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, 4.0])
        node.SetValue(KratosMultiphysics.REACTION, [0.0, 0.0, 4.0])

        node = model_part.CreateNewNode(6, 2.0, 1.0, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [5.0, 0.0, 0.0])
        node.SetValue(KratosMultiphysics.DISPLACEMENT, [5.0, 0.0, 0.0])
        node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, 5.0])
        node.SetValue(KratosMultiphysics.REACTION, [0.0, 0.0, 5.0])

        model_part.CreateNewElement(
            "Element2D4N",
            1,
            [1, 2, 3, 4],
            model_part.GetProperties()[1])

        model_part.CreateNewElement(
            "Element2D4N",
            2,
            [2, 5, 6, 3],
            model_part.GetProperties()[1])

        return model, model_part


    def test_Vertex(self):
        model, model_part = self.MakeModelPart()
        locator = HDF5Application.BruteForcePointLocatorAdaptor(
            model_part,
            KratosMultiphysics.Configuration.Initial,
            1e-6)

        vertices: 'list[HDF5Application.Vertex]' = []
        vertices.append(HDF5Application.Vertex([0.5, 0.5, 0.0], locator, 1))
        vertices.append(HDF5Application.Vertex([1.5, 0.5, 0.0], locator, 2))

        for vertex in vertices:
            self.assertTrue(vertex.IsLocated())

        self.assertVectorAlmostEqual(
            vertices[0].GetValue(KratosMultiphysics.DISPLACEMENT),
            [1.5, 0.0, 0.0])
        self.assertVectorAlmostEqual(
            vertices[0].GetValue(KratosMultiphysics.REACTION),
            [0.0, 0.0, 1.5])

        self.assertVectorAlmostEqual(
            vertices[1].GetValue(KratosMultiphysics.DISPLACEMENT),
            [3.0, 0.0, 0.0])
        self.assertVectorAlmostEqual(
            vertices[1].GetValue(KratosMultiphysics.REACTION),
            [0.0, 0.0, 3.0])

        self.assertVectorAlmostEqual(
            vertices[0].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT),
            [1.5, 0.0, 0.0])
        self.assertVectorAlmostEqual(
            vertices[0].GetSolutionStepValue(KratosMultiphysics.REACTION),
            [0.0, 0.0, 1.5])

        self.assertVectorAlmostEqual(
            vertices[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT),
            [3.0, 0.0, 0.0])
        self.assertVectorAlmostEqual(
            vertices[1].GetSolutionStepValue(KratosMultiphysics.REACTION),
            [0.0, 0.0, 3.0])


if __name__ == "__main__":
    KratosUnittest.main()