import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication as DEM

class TestSearchNodes(KratosUnittest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.base_model_part = self.model.CreateModelPart('Base')
        self.target_model_part = self.model.CreateModelPart('Target')
        # creating a search tool with a periodic bounding box
        dx = self.unit_length = 0.1
        self.epsilon = 0.0001
        self.search_strategy = DEM.OMP_DEMSearch(0, 0, 0, 10*dx, 10*dx, 10*dx)
        self.CreateNodes()

    def test_SearchNodesInTargetModelPart(self):
        dx = self.unit_length
        epsilon = self.epsilon

        neighbours = []
        # setting search radius just about large enough to catch some of the neighbours
        radii = [2*dx + epsilon for node in self.base_model_part.Nodes]
        self.search_strategy.SearchNodesInRadiusExclusive(self.base_model_part.Nodes,
                                                          self.target_model_part.Nodes,
                                                          radii,
                                                          DEM.VectorResultNodesContainer(),
                                                          DEM.VectorDistances(),
                                                          neighbours)

        self.assertEqual(set(neighbours[0]), {3, 4})
        self.assertEqual(set(neighbours[1]), {4})

        neighbours = []
        # setting search radius a tad too short to catch some of the neighbours
        radii = [2*dx - epsilon for node in self.base_model_part.Nodes]
        self.search_strategy.SearchNodesInRadiusExclusive(self.base_model_part.Nodes,
                                                          self.target_model_part.Nodes,
                                                          radii,
                                                          DEM.VectorResultNodesContainer(),
                                                          DEM.VectorDistances(),
                                                          neighbours)

        self.assertEqual(set(neighbours[0]), set())
        self.assertEqual(set(neighbours[1]), {4})

    def CreateNodes(self):
        dx = self.unit_length
        # Base nodes, for which neighbours must be searched for
        # ------------------------------------------
        self.CreateNode('Base', 1, [dx, dx, dx])
        self.CreateNode('Base', 2, [0.9*dx, 0.9*dx, 0.9*dx])

        # Target nodes, which are search candidates
        # ------------------------------------------
        # the first one is simply displaced two length units in the x-direction with respect to
        # the first base node
        self.CreateNode('Target', 3, [3*dx, dx, dx])
        # the second candidate is displaced along the diagonal two length units towards the origin,
        # so that it comes in through the opposite end of the diagonal
        dx_opp = (11 - 2.0/3**0.5) * dx
        self.CreateNode('Target', 4, [dx_opp, dx_opp, dx_opp])
        # the third one is displaced 3 units in the y direction from the first base node
        self.CreateNode('Target', 5, [4*dx, dx, dx])
        # the fourth candidate is also displaced 3 units in the y direction, but in the opposite
        # direction (through the boundary)
        self.CreateNode('Target', 6, [1 - 2*dx, dx, dx])

    def CreateNode(self, model_part_name, id, coordinates):
        model_part = self.model.GetModelPart(model_part_name)
        node = model_part.CreateNewNode(id, coordinates[0], coordinates[1], coordinates[2])


if __name__ == '__main__':
    KratosUnittest.main()
