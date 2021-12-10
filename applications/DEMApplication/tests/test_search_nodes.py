import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication as DEM
import numpy as np

class TestSearchNodes(KratosUnittest.TestCase):

    def setUp(self):
        self.current_model = Kratos.Model()
        self.main_model_part = self.current_model.CreateModelPart("Main")
        self.target_model_part = self.current_model.CreateModelPart("Target")

        # recycling 'ACTIVATION_LEVEL' variable to identify the particle of which
        # the node is a neighbour.
        self.main_model_part.AddNodalSolutionStepVariable(Kratos.RADIUS)
        self.main_model_part.AddNodalSolutionStepVariable(Kratos.ACTIVATION_LEVEL)
        self.target_model_part.AddNodalSolutionStepVariable(Kratos.RADIUS)
        self.target_model_part.AddNodalSolutionStepVariable(Kratos.ACTIVATION_LEVEL)
        # creating a search tool with a periodic bounding box
        self.search_strategy = DEM.OMP_DEMSearch(0, 0, 0, 1, 1, 1)
        self.search_radius = 0.21
        self.CreateNodes()

    def test_SearchNodesInTargetModelPart(self):
        radii = [self.search_radius]
        results = DEM.VectorResultNodesContainer()
        distances = DEM.VectorDistances()
        neighbours = []
        self.search_strategy.SearchNodesInRadiusExclusive(self.main_model_part.Nodes,
                                                          self.target_model_part.Nodes,
                                                          radii,
                                                          results,
                                                          distances,
                                                          neighbours)

        print('neighbours', neighbours)
        # self.assertEqual(counter, 1)

    def CreateNodes(self):
        node1 = self.main_model_part.CreateNewNode(1, 0.1, 0.1, 0.1)
        node2 = self.target_model_part.CreateNewNode(2, 0.2, 0.1, 0.1)
        coors_diagonal = 1 - 0.2/3**0.5
        node3 = self.target_model_part.CreateNewNode(3, coors_diagonal, coors_diagonal, coors_diagonal)
        node4 = self.target_model_part.CreateNewNode(4, 0.2, 0.1, 0.1)

        node1.SetSolutionStepValue(Kratos.RADIUS, self.search_radius)
        node2.SetSolutionStepValue(Kratos.RADIUS, self.search_radius)
        node3.SetSolutionStepValue(Kratos.RADIUS, self.search_radius)
        node4.SetSolutionStepValue(Kratos.RADIUS, self.search_radius)


        counter = 0
        for node in self.target_model_part.Nodes:
            counter += 1

        self.assertEqual(counter, 3)


    def test_CreateWallQuadrilateral(self):
        pass


if __name__ == '__main__':
    KratosUnittest.main()
